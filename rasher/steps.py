"""Per-sample steps: bait, filter and assemble reads."""

from __future__ import annotations

import logging
import re
import shutil
import tempfile
from contextlib import nullcontext
from dataclasses import dataclass, field
from itertools import zip_longest
from pathlib import Path

from rasher import RasherError
from rasher.samples import Sample
from rasher.seqio import (
    ReadStats,
    Record,
    acgtn,
    canonical_name,
    concatenate,
    extract_pairs,
    extract_reads,
    n50,
    open_write,
    read_records,
    read_stats,
    split_extension,
    subsample_reads,
    write_fasta,
)
from rasher.tools import ToolError, run, which

log = logging.getLogger(__name__)


class SampleFailed(RasherError):
    """This sample cannot go further; the others can. `stats`: what is known about it (e.g. its raw reads)."""

    def __init__(self, message: str, stats: dict[str, object] | None = None):
        super().__init__(message)
        self.stats = stats or {}


@dataclass
class StepResult:
    output: Path | None
    output2: Path | None = None  # The R2 reads of a paired output (bait, filter); None otherwise
    stats: dict[str, object] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def reads_suffix(fmt: str) -> str:
    return ".fastq.gz" if fmt == "fastq" else ".fasta.gz"


def _read_outputs(out_dir: Path, name: str, fmt: str, paired: bool) -> list[Path]:
    """<sample>_R1 and <sample>_R2 reads files when paired, else <sample>."""
    return [out_dir / f"{name}{tag}{reads_suffix(fmt)}" for tag in (("_R1", "_R2") if paired else ("",))]


# ---------------------------------------------------------------------------------------------------------------
# 1. Baiting
# ---------------------------------------------------------------------------------------------------------------

def _names_from_paf(paf: Path) -> set[str]:
    """The canonical names (without /1 or /2) of the reads aligned in a PAF file."""
    with open(paf, encoding="ascii", errors="replace") as fh:  # As the reads are read (seqio)
        return {canonical_name(line.split("\t", 1)[0]) for line in fh if line.strip()}


def bait_minimap2(sample: Sample, reference: Path, out_dir: Path, log_dir: Path, threads: int,
                  keep_bam: bool) -> StepResult:
    """Keep the reads with at least one alignment to the reference; of paired reads, the pairs with a mate that
    has one.

    Each file is aligned on its own, as single reads: minimap2's short-read mode takes at most two query files,
    and pairs them (which does not change whether a read aligns). With `keep_bam`, the baited reads (not all the
    reads) are aligned again, as pairs when paired, into <sample>.bam.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample.name}.log"
    options = ["sr", "--secondary=no", "-t", str(threads), str(reference)]
    outputs = _read_outputs(out_dir, sample.name, sample.fmt, bool(sample.pairs))
    bam, paf = out_dir / f"{sample.name}.bam", out_dir / f".{sample.name}.paf"
    try:
        names: set[str] = set()
        for path in sample.files:
            run(["minimap2", "-x", *options, str(path)], log_file, stdout=paf, what=f"(mapping {sample.name})")
            names |= _names_from_paf(paf)
        if sample.pairs:
            total, kept = extract_pairs(sample.pairs, names, outputs[0], outputs[1])
        else:
            total, kept = extract_reads(sample.files, names, outputs[0])
        result = _bait_result(outputs, total, kept)
        if keep_bam:
            run([["minimap2", "-ax", *options, *map(str, outputs)],
                 ["samtools", "sort", "-@", str(max(1, threads // 4)), "-o", str(bam), "-"]], log_file,
                what=f"(mapping {sample.name})")
            run(["samtools", "index", str(bam)], log_file)
        return result
    except BaseException:
        # A failed (or interrupted) sample leaves no reads, partial or stale
        for path in (*outputs, *((bam, out_dir / f"{bam.name}.bai") if keep_bam else ())):
            path.unlink(missing_ok=True)
        raise
    finally:
        paf.unlink(missing_ok=True)


def _extension(path: Path) -> str:
    """The extension of a sequence file, compression included ('.fastq.gz'), as its name gives it."""
    parts = split_extension(path.name)
    return path.name[len(parts[0]):] if parts else "".join(path.suffixes)


def _content_extension(path: Path, fmt: str) -> str:
    """'.fastq.gz', '.fasta', ...: the extension telling BBTools (which goes by the extension, case-sensitively)
    what a reads file holds, from its content (as RASHER reads it), whatever its own name says."""
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    return f".{fmt}{'.gz' if gzipped else ''}"


_BBDUK_COUNTS = re.compile(r"^(Input|Contaminants):\s+(\d+) reads\s.*?(\d+) bases", re.M)


def bait_bbduk(sample: Sample, reference: Path, out_dir: Path, log_dir: Path, threads: int, kmer: int,
               hdist: int, memory_gb: int) -> StepResult:
    """Keep the reads sharing at least one k-mer (up to `hdist` mismatches) with the reference; of paired reads,
    the pairs with a mate that does (BBDuk's default for pairs)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample.name}.log"
    paired = bool(sample.pairs)
    tags = ("_R1", "_R2") if paired else ("",)
    mates = [[pair[i] for pair in sample.pairs] for i in (0, 1)] if paired else [sample.files]
    outputs = _read_outputs(out_dir, sample.name, sample.fmt, paired)
    # BBDuk takes one input file (of each mate): several are merged first
    merged = [out_dir / f".{sample.name}.merged{tag}{reads_suffix(sample.fmt)}" if len(files) > 1 else None
              for tag, files in zip(tags, mates)]
    start = log_file.stat().st_size if log_file.exists() else 0
    # bbduk.sh passes its arguments through `eval`: a path with a space (or a shell character) breaks it. BBDuk
    # runs in out_dir, on links named after the sample.
    links = {"ref": out_dir / f".{sample.name}.reference.fasta"}
    targets = {"ref": reference}
    try:
        for key, tag, files, merge in zip(("in", "in2"), tags, mates, merged):
            if merge is not None:
                concatenate(files, merge)
            targets[key] = merge or files[0]
            links[key] = out_dir / f".{sample.name}.input{tag}{_content_extension(targets[key], sample.fmt)}"
        for key, link in links.items():
            link.unlink(missing_ok=True)
            link.symlink_to(targets[key].resolve())
        run(["bbduk.sh", f"-Xmx{max(1, memory_gb)}g", "-eoom", "overwrite=true",
             *(f"{key}={links[key].name}" for key in ("in", "in2") if key in links), f"ref={links['ref'].name}",
             f"threads={threads}", f"k={kmer}", f"hdist={hdist}", "maskmiddle=f",
             *(f"{key}={out.name}" for key, out in zip(("outm", "outm2"), outputs))], log_file, cwd=out_dir,
            what=f"(baiting {sample.name})")
        with open(log_file) as fh:
            fh.seek(start)
            counts = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in _BBDUK_COUNTS.finditer(fh.read())}
        total = ReadStats(*counts["Input"]) if "Input" in counts else None
        return _bait_result(outputs, total, read_stats(outputs))
    except BaseException as exc:
        for out in outputs:  # A failed (or interrupted) sample leaves no reads, partial or stale
            out.unlink(missing_ok=True)
        if isinstance(exc, ToolError):
            with open(log_file, errors="replace") as fh:
                fh.seek(start)
                # Not "OutOfMemoryError" alone: bbduk.sh first prints its java command, with
                # -XX:+ExitOnOutOfMemoryError
                if "java.lang.OutOfMemoryError" in fh.read():
                    raise SampleFailed(f"BBDuk ran out of memory with {memory_gb} GB: increase -m, lower -p, or "
                                       f"use --hdist {max(0, hdist - 1)} (each mismatch allowed multiplies the "
                                       f"memory needed); see {log_file}") from exc
        raise
    finally:
        for link in links.values():
            link.unlink(missing_ok=True)
        for merge in merged:
            if merge is not None:
                merge.unlink(missing_ok=True)


def _bait_result(outputs: list[Path], total: ReadStats | None, kept: ReadStats) -> StepResult:
    """`outputs`: the reads file, or the R1 and R2 files. The statistics count the reads of both mates."""
    stats: dict[str, object] = {"Raw_reads": total.reads if total else "NA",
                                "Raw_bases": total.bases if total else "NA",
                                "Baited_reads": kept.reads, "Baited_bases": kept.bases}
    stats["Baited_pct"] = (f"{100 * kept.bases / total.bases:.3f}" if total and total.bases else "NA")
    if kept.reads == 0:
        raise SampleFailed("no reads matched the reference", stats)
    return StepResult(outputs[0], outputs[1] if len(outputs) > 1 else None, stats)


# ---------------------------------------------------------------------------------------------------------------
# 2. Read filtering
# ---------------------------------------------------------------------------------------------------------------

def filter_fastp(name: str, reads: Path, reads2: Path | None, fmt: str, out_dir: Path, log_dir: Path, *,
                 min_length: int, target_bases: int, threads: int) -> StepResult:
    """Trim and filter the reads with fastp (its defaults: adapters, base qualities; a pair goes when either
    mate does), then keep an evenly spread subset of about `target_bases` (`subsample_reads`).

    fastp cannot read fasta: fasta reads (or pairs) shorter than `min_length` are dropped here instead, without
    trimming. fastp's report (<sample>.fastp.json, for MultiQC, and .html) goes to `log_dir`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}.log"
    paired = reads2 is not None
    outputs = _read_outputs(out_dir, name, fmt, paired)
    trimmed = _read_outputs(out_dir, f".{name}.trimmed", fmt, paired)
    out2, trimmed2 = (outputs[1], trimmed[1]) if paired else (None, None)
    try:
        if fmt == "fasta":
            _drop_short(reads, reads2, trimmed[0], trimmed2, min_length)
        else:
            mates = ["--in2", str(reads2), "--out2", str(trimmed2)] if paired else []
            run(["fastp", "--in1", str(reads), "--out1", str(trimmed[0]), *mates, "--length_required",
                 str(min_length), "--thread", str(threads), "--json", str(log_dir / f"{name}.fastp.json"),
                 "--html", str(log_dir / f"{name}.fastp.html")], log_file, what=f"(filtering {name})")
        subsample_reads(trimmed[0], trimmed2, outputs[0], out2, target_bases)
        stats = read_stats(outputs)
        if stats.reads == 0:
            raise SampleFailed(f"no reads left after filtering (all shorter than {min_length} bp"
                               f"{'' if fmt == 'fasta' else ' or of low quality'}?)")
        return StepResult(outputs[0], out2, {"Filtered_reads": stats.reads, "Filtered_bases": stats.bases,
                                             "Filtered_N50": stats.n50})
    except BaseException:
        for out in outputs:  # A failed (or interrupted) sample leaves no reads, partial or stale
            out.unlink(missing_ok=True)
        raise
    finally:
        for tmp in trimmed:
            tmp.unlink(missing_ok=True)


def _drop_short(reads: Path, reads2: Path | None, out: Path, out2: Path | None, min_length: int) -> None:
    """Copy the reads of at least `min_length` bp; of pairs (`reads2`), those whose mates both are."""
    first, second = read_records(reads), (read_records(reads2) if reads2 is not None else None)
    with open_write(out) as o1, (open_write(out2) if out2 is not None else nullcontext()) as o2:
        for rec1, rec2 in zip_longest(first, second or ()):
            if rec1 is None or (second is not None and rec2 is None):
                raise RasherError(f"{reads} and {reads2} do not hold the same number of reads")
            if len(rec1.seq) >= min_length and (rec2 is None or len(rec2.seq) >= min_length):
                o1.write(rec1.format())
                if rec2 is not None:
                    o2.write(rec2.format())


# ---------------------------------------------------------------------------------------------------------------
# 3. Assembly
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class AssemblyDirs:
    root: Path  # 3_assembled

    @property
    def assemblies(self) -> Path:
        return self.root / "all_assemblies"

    @property
    def graphs(self) -> Path:
        return self.root / "assembly_graphs"

    def work(self, name: str) -> Path:
        return self.root / name

    def clear_work(self, name: str) -> Path:
        """The sample's working folder, deleted (the files of an earlier assembly). Refused unless the folder is
        right in `root` (sample names are checked, but a name such as '..' would delete the output folder)."""
        work = self.work(name)
        if work.resolve().parent != self.root.resolve():
            raise SampleFailed(f"refusing to delete {work}: not a sample folder in {self.root}")
        shutil.rmtree(work, ignore_errors=True)
        return work

    def clear_graph(self, name: str) -> None:
        """The graph of an earlier assembly of the sample, deleted (this assembler may draw none)."""
        for old in (self.graphs / f"{name}.gfa", self.graphs / f"{name}.png"):
            old.unlink(missing_ok=True)


def _publish_assembly(name: str, contigs: Path, dirs: AssemblyDirs, *, rename: bool,
                      circular: set[str] | None = None) -> list[Record]:
    """Copy the contigs to all_assemblies/<sample>.fasta, prefixing the contig names with the sample name.

    Contigs known to be circular (in `circular`, or flagged so in their header by the assembler) are marked
    `circular=true` in their new header.
    """
    records = [r for r in read_records(contigs) if r.seq]
    if not records:
        raise SampleFailed("the assembler produced no contigs")
    dirs.assemblies.mkdir(parents=True, exist_ok=True)
    out = []
    for i, rec in enumerate(records, 1):
        header = f"{name}_{i}" if rename else f"{name}_{rec.name}"
        if rec.name in (circular or set()) or _CIRCULAR_TAGS.search(rec.header):
            header += " circular=true"
        out.append(Record(header, acgtn(rec.seq)))  # Ambiguity codes (e.g. from consensus calls) become N
    write_fasta(dirs.assemblies / f"{name}.fasta", out)
    return out


def _circular_count(records: list[Record]) -> int:
    return sum(1 for r in records if r.header.endswith(" circular=true"))


_CIRCULAR_TAGS = re.compile(r"(?:\bXO:i:1\b|circular[-=](?:yes|true|y)(?![a-z]))", re.I)


def count_circular(contigs: Path) -> int:
    """Contigs flagged circular in their header (myloasm: circular-yes; Raven: XO:i:1)."""
    return sum(1 for rec in read_records(contigs) if _CIRCULAR_TAGS.search(rec.header))


def _publish_graph(name: str, gfa: Path, dirs: AssemblyDirs, log_file: Path) -> None:
    if not gfa.is_file():
        return
    dirs.graphs.mkdir(parents=True, exist_ok=True)
    target = dirs.graphs / f"{name}.gfa"
    shutil.copyfile(gfa, target)
    if which("Bandage"):
        try:  # A picture of the graph is nice to have; its failure does not fail the sample
            run(["Bandage", "image", str(target), str(target.with_suffix(".png"))], log_file,
                env={"QT_QPA_PLATFORM": "offscreen"})
        except RasherError as exc:
            log.warning("Bandage could not draw the graph of %s: %s", name, str(exc).splitlines()[0])


def _assembly_stats(records: list[Record], reference_length: int) -> dict[str, object]:
    lengths = [len(r.seq) for r in records]
    total = sum(lengths)
    return {"Contigs": len(lengths), "Assembly_length": total, "Largest_contig": max(lengths),
            "Assembly_N50": n50(lengths), "Length_vs_reference": f"{total / reference_length:.3f}"}


def _denovo_stats(records: list[Record], reference_length: int) -> dict[str, object]:
    """Assembly statistics of SPAdes and SKESA: both give k-mer coverage, not read depth (NA)."""
    return {**_assembly_stats(records, reference_length), "Circular_contigs": _circular_count(records),
            "Assembly_depth": "NA"}


_SHELL_CHARACTERS = re.compile(r"""[\s'"`$\\!*?\[\](){}<>|&;#,]""")


def _shell_safe(*paths: Path) -> bool:
    """Whether the absolute paths have no whitespace or shell character, nor a comma: what a program that writes
    paths into configuration files, runs commands through a shell, or splits lists on commas may break on."""
    return not any(_SHELL_CHARACTERS.search(str(p.absolute())) for p in paths)


def assemble_spades(name: str, reads: Path, reads2: Path | None, dirs: AssemblyDirs, log_dir: Path, *,
                    threads: int, memory_gb: int, reference_length: int) -> StepResult:
    """De novo assembly with SPAdes --isolate (which skips read correction: fasta reads assemble as well). SPAdes
    flags no contig circular, and the _cov_ of its contig names is k-mer coverage: no Assembly_depth."""
    work = dirs.clear_work(name)  # SPAdes wants a new (or empty) output folder
    dirs.clear_graph(name)
    dirs.root.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}.log"
    mates = [("-1", reads), ("-2", reads2)] if reads2 is not None else [("-s", reads)]
    # SPAdes writes its paths into configuration files that a space breaks: with such paths, it reads and writes
    # through links in a temporary folder
    link_dir = (Path(tempfile.mkdtemp(prefix="rasher_spades_"))
                if not _shell_safe(work, *(path for _, path in mates)) else None)
    try:
        out = work
        if link_dir is not None:
            if not _shell_safe(link_dir):
                raise SampleFailed(f"SPAdes cannot use paths with spaces or shell characters ({work}), and the "
                                   f"temporary folder {link_dir.parent} has some too: set TMPDIR to another "
                                   "folder")
            work.mkdir()
            out = link_dir / "out"
            out.symlink_to(work.resolve())
            links = [(option, link_dir / f"reads{option}{_extension(path)}") for option, path in mates]
            for (_, link), (_, path) in zip(links, mates):
                link.symlink_to(path.resolve())
            mates = links
        run(["spades.py", "--isolate", *(arg for option, path in mates for arg in (option, str(path))),
             "-t", str(threads), "-m", str(memory_gb), "-o", str(out)], log_file, what=f"(assembling {name})")
    finally:
        if link_dir is not None:
            shutil.rmtree(link_dir, ignore_errors=True)  # The links only
    contigs = work / "contigs.fasta"
    if not contigs.is_file():
        raise SampleFailed(f"SPAdes produced no assembly (see {work / 'spades.log'})")
    records = _publish_assembly(name, contigs, dirs, rename=False)  # NODE_1_length_..._cov_...: long, unique
    _publish_graph(name, work / "assembly_graph_with_scaffolds.gfa", dirs, log_file)
    return StepResult(dirs.assemblies / f"{name}.fasta", stats=_denovo_stats(records, reference_length))


def assemble_skesa(name: str, reads: Path, reads2: Path | None, dirs: AssemblyDirs, log_dir: Path, *,
                   threads: int, memory_gb: int, reference_length: int) -> StepResult:
    """De novo assembly with SKESA. Conservative, SKESA ends contigs at repeats rather than guess through them:
    an organelle's inverted repeats split its assembly into several contigs. No graph; no Assembly_depth (its
    contig names give k-mer coverage)."""
    work = dirs.clear_work(name)
    dirs.clear_graph(name)
    work.mkdir(parents=True)
    log_file = log_dir / f"{name}.log"
    # A comma in --reads separates the mates: SKESA reads links in its folder, and runs from there
    mates = [("", reads)] if reads2 is None else [("_R1", reads), ("_R2", reads2)]
    links = {work / f"reads{tag}{_extension(path)}": path for tag, path in mates}
    try:
        for link, path in links.items():
            link.symlink_to(path.resolve())
        run(["skesa", "--reads", ",".join(link.name for link in links), "--contigs_out", "contigs.fasta",
             "--cores", str(threads), "--memory", str(memory_gb)], log_file, cwd=work,
            what=f"(assembling {name})")
    finally:
        for link in links:
            link.unlink(missing_ok=True)
    contigs = work / "contigs.fasta"
    if not contigs.is_file() or not any(rec.seq for rec in read_records(contigs)):
        raise SampleFailed(f"SKESA produced no contigs (see {log_file})")
    records = _publish_assembly(name, contigs, dirs, rename=False)  # Contig_1_42.5: number, k-mer coverage
    return StepResult(dirs.assemblies / f"{name}.fasta", stats=_denovo_stats(records, reference_length))


def off_reference(contigs: Path, reference: Path, log_file: Path, threads: int) -> float:
    """Percentage of a de novo assembly's bases that align nowhere on the reference (minimap2 -x asm20): contigs
    of another origin (nuclear DNA, another organelle) or sequence the reference lacks."""
    lengths = {r.name: len(r.seq) for r in read_records(contigs)}
    total = sum(lengths.values())
    if not total:
        return 0.0
    with tempfile.TemporaryDirectory(prefix="rasher_offref_") as tmp:
        paf = Path(tmp) / "contigs.paf"
        run(["minimap2", "-x", "asm20", "-t", str(threads), str(reference), str(contigs)], log_file, stdout=paf,
            what="(aligning the contigs to the reference)")
        spans: dict[str, list[tuple[int, int]]] = {}
        with open(paf) as fh:
            for line in fh:
                fields = line.split("\t")
                if len(fields) >= 12 and fields[0] in lengths:
                    spans.setdefault(fields[0], []).append((int(fields[2]), int(fields[3])))
    aligned = 0
    for intervals in spans.values():
        end = 0
        for start, stop in sorted(intervals):
            aligned += max(0, stop - max(start, end))
            end = max(end, stop)
    return 100 * (total - aligned) / total


def assemble_samtools(name: str, reads: Path, reads2: Path | None, reference: Path, dirs: AssemblyDirs,
                      log_dir: Path, *, threads: int, reference_length: int, fill_gaps: bool,
                      cutoff: int = 20) -> StepResult:
    """Templated assembly: samtools consensus of the reads aligned to the reference with minimap2 (as pairs when
    `reads2` is given: two query files are mates to minimap2)."""
    work = dirs.clear_work(name)
    dirs.clear_graph(name)  # From an earlier de novo assembly
    work.mkdir(parents=True)
    log_file = log_dir / f"{name}.log"
    bam = work / "aligned.bam"
    queries = [str(reads)] + ([str(reads2)] if reads2 is not None else [])
    run([["minimap2", "-ax", "sr", "--secondary=no", "-t", str(threads), str(reference), *queries],
         ["samtools", "sort", "-@", str(max(1, threads // 4)), "-o", str(bam), "-"]], log_file,
        what=f"(aligning {name})")
    run(["samtools", "index", str(bam)], log_file)
    consensus = work / "consensus.fasta"
    # -a: every position of the reference sequences that have reads, N where no read (a reference sequence
    # without any read is left out); Illumina error profile; insertions shown, deletions removed, so the
    # consensus follows the sample, not the reference; -C: a call below this quality is N.
    run(["samtools", "consensus", "-a", "-X", "hiseq", "--show-ins", "yes", "--show-del", "no",
         "-d", "3", "-C", str(cutoff), "-f", "fasta", "-o", str(consensus), str(bam)], log_file,
        what=f"(consensus {name})")
    records = [r for r in read_records(consensus) if r.seq]
    if fill_gaps:  # Replace N runs by the reference: only safe when the sample is known to match it there
        refs = {r.name: r.seq for r in read_records(reference)}
        records = [Record(r.header, _fill_from(r.seq, refs.get(r.name, ""))) for r in records]
        write_fasta(consensus, records)
    published = _publish_assembly(name, consensus, dirs, rename=False)
    stats = _assembly_stats(published, reference_length)
    n_bases = sum(r.seq.count("N") for r in published)
    stats["N_bases"] = n_bases
    notes = [f"{n_bases:,} N bases (no or ambiguous read support)"] if n_bases else []
    missing = [r.name for r in read_records(reference) if r.name not in {c.name for c in records}]
    if missing:
        notes.append(f"no reads on reference sequence(s) {', '.join(missing)}: not in the assembly")
    return StepResult(dirs.assemblies / f"{name}.fasta", stats=stats, notes=notes)


def _fill_from(seq: str, ref: str) -> str:
    """Leading/trailing N runs map to the reference ends; inner ones cannot be placed reliably, keep them."""
    if not ref:
        return seq
    core = seq.strip("N")
    if not core:  # No read support at all: nothing to anchor the reference on
        return seq
    lead = len(seq) - len(seq.lstrip("N"))
    trail = len(seq) - len(seq.rstrip("N"))
    return ref[:lead] + core + (ref[len(ref) - trail:] if trail else "")
