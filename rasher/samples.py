"""Finding the samples to process: an input file, an input folder, or a sample sheet."""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from rasher import RasherError
from rasher.metadata import read_table
from rasher.seqio import sniff_format, split_extension

log = logging.getLogger(__name__)

# Folders that MinKNOW/Dorado write next to the barcode folders and that are not samples.
SKIPPED_FOLDERS = {"unclassified", "mixed"}
# The whole name (fullmatch). Not starting with a dot: "." and ".." would name the step folder or its parent (a
# sample's assembly folder is deleted before each assembly), and the files of hidden names are skipped in folders.
VALID_NAME = re.compile(r"[A-Za-z0-9_+-][A-Za-z0-9._+-]*")
NAME_RULE = "use only letters, digits and . _ + -, not starting with a dot"


@dataclass
class Sample:
    name: str
    files: list[Path]  # All the files used: when paired, the pairs flattened (R1a, R2a, R1b, R2b, ...)
    fmt: str = "fastq"  # 'fastq' or 'fasta'
    pairs: list[tuple[Path, Path]] = field(default_factory=list)  # (R1, R2) pairs; empty = single-end


# The mate markers before the extension: X_R1/X_R2, X_1/X_2, X.R1/X.R2 (R or r), with bcl2fastq's _001 tail
# allowed (X_S1_L001_R1_001.fastq.gz). The base is greedy: of "a_1_2", the final _2 is the marker.
PAIR_STEM = re.compile(r"(?P<base>.+)(?P<sep>[._])(?P<r>[Rr]?)(?P<mate>[12])(?P<tail>_001)?$")


def split_pair_stem(filename: str) -> tuple[str, tuple[str, str, str], int] | None:
    """(base, pairing key, mate 1 or 2) when the file name ends in a mate marker before its fasta/fastq
    extension, else None. Mates pair when their base and key match ('R1' pairs with 'r2'; '_1' does not pair
    with '_R2')."""
    parts = split_extension(filename)
    if not parts:
        return None
    match = PAIR_STEM.fullmatch(parts[0])
    if not match:
        return None
    return match["base"], (match["sep"], match["r"].upper(), match["tail"] or ""), int(match["mate"])


def pair_files(files: list[Path]) -> tuple[list[tuple[Path, Path]], list[Path]]:
    """The (R1, R2) pairs among `files`, and the files left unpaired (in their original order).

    Two files pair when they sit in the same folder and their names differ only by the mate marker (the format
    their extension gives included: X_R1.fasta does not pair with X_R2.fastq), and no third file claims the same
    marker: two X_R1 files leave all of X unpaired (a duplicate to report, not to guess around). A file and a
    link to it do not pair."""
    groups: dict[tuple[object, ...], dict[int, list[Path]]] = {}
    for f in files:
        if split := split_pair_stem(f.name):
            base, key, mate = split
            fmt = split_extension(f.name)[1]  # type: ignore[index] - a pair stem has an extension
            groups.setdefault((f.parent.resolve(), base, key, fmt), {}).setdefault(mate, []).append(f)
    pairs, paired = [], set()
    for mates in groups.values():
        if len(mates.get(1, [])) == 1 and len(mates.get(2, [])) == 1:
            r1, r2 = mates[1][0], mates[2][0]
            if r1.resolve() == r2.resolve():
                log.warning("%s and %s are the same file (one links to the other): not paired", r1, r2)
                continue
            pairs.append((r1, r2))
            paired.update((r1, r2))
    pairs.sort()
    return pairs, [f for f in files if f not in paired]


def _sequence_files(folder: Path) -> list[Path]:
    """The sequence files in a folder and its subfolders, hidden ones (and those in hidden folders) excepted.
    Linked subfolders are followed (Path.rglob does not follow them before Python 3.13), each folder once: a
    link to a folder already read (or to a folder above it) is skipped."""
    files, seen = [], set()
    for root, folders, names in os.walk(folder, followlinks=True):
        real = os.path.realpath(root)
        if real in seen:
            folders[:] = []
            continue
        seen.add(real)
        folders[:] = sorted(f for f in folders if not f.startswith("."))  # In order: a folder's first link wins
        files += [path for name in names if not name.startswith(".") and split_extension(name)
                  and (path := Path(root) / name).is_file()]
    return sorted(files)


def _check_name(name: str, origin: object) -> None:
    if not VALID_NAME.fullmatch(name):
        raise RasherError(f"Invalid sample name {name!r} (from {origin}): {NAME_RULE}")
    if name.lower() in {"reference", "all_assemblies", "assembly_graphs"}:
        raise RasherError(f"Sample name {name!r} (from {origin}) is reserved; rename the file or folder")


def _make_sample(name: str, files: list[Path], origin: object) -> Sample | None:
    _check_name(name, origin)
    unique: dict[Path, Path] = {}
    for f in files:
        if (real := f.resolve()) in unique:
            log.warning("Sample %s: %s is the file %s (a link): its reads are used once", name, f, unique[real])
        else:
            unique[real] = f
    files = list(unique.values())
    formats = {f: sniff_format(f) for f in files}
    empty = [f for f, fmt in formats.items() if fmt is None]
    for f in empty:
        log.warning("Skipping empty file %s", f)
    files = [f for f in files if formats[f] is not None]
    if not files:
        log.warning("Skipping sample %s: no reads", name)
        return None
    kinds = {formats[f] for f in files}
    if len(kinds) > 1:
        raise RasherError(f"Sample {name} mixes fasta and fastq files: {', '.join(map(str, files))}")
    pairs, unpaired = pair_files(files)
    if pairs and unpaired:
        # Pairs and leftovers do not mix: every step would need both a paired and an unpaired stream
        log.warning("Sample %s: %s %s no mate next to the %d pair(s): skipped (list it as its own sample to "
                    "use it)", name, ", ".join(f.name for f in unpaired),
                    "has" if len(unpaired) == 1 else "have", len(pairs))
    if pairs:
        files = [f for pair in pairs for f in pair]
    else:
        marked = [f for f in files if split_pair_stem(f.name)]
        if {split_pair_stem(f.name)[2] for f in marked} == {1, 2}:  # type: ignore[index] - marked
            log.warning("Sample %s: %s look like read mates but were not paired (a duplicate mate file, or mates "
                        "in different folders): processed as single-end reads", name,
                        ", ".join(f.name for f in marked))
    return Sample(name, files, kinds.pop(), pairs)


def discover(path: Path) -> list[Sample]:
    """Samples from a single file, or from a folder.

    In a folder, two files that differ only by a mate marker (X_R1/X_R2, X_1/X_2, X.R1/X.R2 before the
    extension) are one paired sample named by what is before the marker; every other sequence file is one
    single-end sample named after the file; each subfolder is one sample named after the subfolder, with all
    the sequence files it contains (pairs found the same way). In a folder, a file named as a pair's sample (X
    next to X_1 and X_2, as fasterq-dump --split-3 writes its unpaired reads) is skipped.
    """
    if not path.exists():
        raise RasherError(f"Input not found: {path}")
    candidates: list[tuple[str, list[Path], object]] = []
    if path.is_file():
        parts = split_extension(path.name)
        if not parts:
            raise RasherError(f"Not a fasta/fastq file (by extension): {path}")
        if split_pair_stem(path.name):
            log.info("%s looks like one mate of a pair: give the folder or a sample sheet to pair it with its "
                     "mate", path)
        candidates.append((parts[0], [path], path))
    else:
        flat: list[Path] = []
        for entry in sorted(path.iterdir()):
            if entry.name.startswith("."):
                continue
            if entry.is_dir():
                files = _sequence_files(entry)
                if not files:
                    log.info("Skipping folder %s: no fasta/fastq files", entry)
                    continue
                if entry.name in SKIPPED_FOLDERS:
                    log.info("Skipping folder %s (reads without a barcode)", entry)
                    continue
                candidates.append((entry.name, files, entry))
            elif entry.is_file() and split_extension(entry.name):
                flat.append(entry)
        pairs, unpaired = pair_files(flat)
        bases: dict[str, str] = {}
        for r1, r2 in pairs:
            base, _, _ = split_pair_stem(r1.name)  # type: ignore[misc] - r1 matched in pair_files
            bases[base] = f"{r1.name} + {r2.name}"
            candidates.append((base, [r1, r2], bases[base]))
        for f in unpaired:
            stem = split_extension(f.name)[0]  # type: ignore[index] - flat files have an extension
            if stem in bases:
                log.warning("Skipping %s: the sample %s is the pair %s, and this file the unpaired reads "
                            "fasterq-dump writes next to a pair; list it as its own sample to use it", f, stem,
                            bases[stem])
                continue
            candidates.append((stem, [f], f))
        candidates.sort(key=lambda c: c[0])
    return _finalize(candidates, path)


def read_sample_sheet(path: Path) -> list[Sample]:
    """Samples from a TSV or CSV file with the columns 'sample' and 'file' (files separated by ';').

    Several rows may share a sample name; their files are merged. Relative paths are resolved from the
    sheet's folder. Lines starting with '#' are comments, before or after the header (a sample can be left out
    by commenting its line). The other columns are metadata (rasher.metadata.sheet_metadata).
    """
    header, rows = read_table(path, "Sample sheet", comment_lines=True)
    columns = {c.lower(): c for c in header}
    if "sample" not in columns or "file" not in columns:
        raise RasherError(f"Sample sheet {path} needs the columns 'sample' and 'file' (found: {header})")
    grouped: dict[str, list[Path]] = {}
    for row in rows:
        name = row.get(columns["sample"], "")
        files = [f.strip() for f in row.get(columns["file"], "").split(";") if f.strip()]
        if not name or not files:
            raise RasherError(f"Sample sheet {path}: row without a sample name or a file: {row}")
        for f in files:
            file = Path(f).expanduser()
            file = file if file.is_absolute() else path.parent / file
            if not file.is_file():
                raise RasherError(f"Sample sheet {path}: file not found for {name}: {file}")
            files_of = grouped.setdefault(name, [])
            if file.resolve() in {f.resolve() for f in files_of}:
                log.warning("Sample sheet %s: %s is listed twice for %s; its reads are used once",
                            path, file, name)
                continue
            files_of.append(file)
    candidates = [(name, files, path) for name, files in grouped.items()]
    return _finalize(candidates, path)


def _finalize(candidates: list[tuple[str, list[Path], object]], source: Path) -> list[Sample]:
    seen: dict[str, object] = {}
    for name, _, origin in candidates:
        if name in seen:
            raise RasherError(f"Two inputs give the same sample name {name!r}: {seen[name]} and {origin}")
        seen[name] = origin
    samples = [s for name, files, origin in candidates if (s := _make_sample(name, files, origin))]
    if not samples:
        raise RasherError(f"No fasta/fastq files with reads found in {source}")
    names = {s.name for s in samples}
    for s in samples:  # A paired sample X writes X_R1/X_R2 files, which a sample of that name would share
        if s.pairs:
            for mate in ("_R1", "_R2"):
                if s.name + mate in names:
                    raise RasherError(f"Sample names {s.name!r} (paired) and {s.name + mate!r} collide: the "
                                      f"paired sample writes {s.name + mate} files; rename one")
        elif s.name + "_R1" in names and s.name + "_R2" in names:
            # A single-end sample deletes the X_R1 + X_R2 files of an earlier paired run of it (rasher.steps)
            raise RasherError(f"Sample names {s.name!r}, {s.name + '_R1'!r} and {s.name + '_R2'!r} collide: the "
                              f"reads files of the last two are those {s.name} writes when paired; rename one")
    return samples
