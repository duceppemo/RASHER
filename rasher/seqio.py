"""Reading and writing fasta/fastq files, gzipped or not."""

from __future__ import annotations

import gzip
import io
import re
import zlib
from collections.abc import Iterable, Iterator
from contextlib import nullcontext
from dataclasses import dataclass, field
from itertools import zip_longest
from pathlib import Path

from rasher import RasherError
from rasher.tools import stopping

# Errors raised while decompressing a truncated or corrupt gzip file.
DECOMPRESSION_ERRORS: tuple[type[BaseException], ...] = (EOFError, zlib.error, gzip.BadGzipFile)
try:  # Optional: python-isal decompresses gzip several times faster
    from isal import igzip as _gzip
    from isal.igzip_lib import IsalError

    DECOMPRESSION_ERRORS += (IsalError,)
except ImportError:  # pragma: no cover - depends on the environment
    _gzip = gzip

FASTQ_EXTENSIONS = (".fastq", ".fq")
FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas")


def split_extension(filename: str) -> tuple[str, str] | None:
    """Return (stem, 'fastq' or 'fasta') for a sequence file name, or None if it is not one."""
    lower = filename.lower()
    for fmt, extensions in (("fastq", FASTQ_EXTENSIONS), ("fasta", FASTA_EXTENSIONS)):
        for ext in extensions:
            for suffix in (ext + ".gz", ext):
                if lower.endswith(suffix) and len(filename) > len(suffix):
                    return filename[: -len(suffix)], fmt
    return None


UTF8_BOM = b"\xef\xbb\xbf"


def open_text(path: Path) -> io.TextIOBase:
    """Open a text file for reading, decompressing it if it is gzipped (detected from its content), without the
    UTF-8 byte order mark that some editors write at its start (it would hide the first line's '>' or '@')."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    binary = _gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")  # noqa: SIM115 - returned
    try:
        if binary.read(len(UTF8_BOM)) != UTF8_BOM:
            binary.seek(0)
    except BaseException:
        binary.close()
        raise
    return io.TextIOWrapper(binary, encoding="ascii", errors="replace")


def open_write(path: Path) -> io.TextIOBase:
    """Open a text file for writing, gzipped if its name ends with .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", compresslevel=4, encoding="ascii", errors="replace")
    return open(path, "w", encoding="ascii", errors="replace")


@dataclass
class Record:
    header: str  # Header line without the leading '>' or '@'
    seq: str
    qual: str | None = None

    @property
    def name(self) -> str:
        return self.header.split(maxsplit=1)[0] if self.header else ""

    def format(self) -> str:
        if self.qual is None:
            return f">{self.header}\n{self.seq}\n"
        return f"@{self.header}\n{self.seq}\n+\n{self.qual}\n"


def sniff_format(path: Path) -> str | None:
    """'fastq', 'fasta', or None for an empty file, from the first non-blank character."""
    try:
        return _sniff_format(path)
    except DECOMPRESSION_ERRORS as exc:
        raise RasherError(f"{path}: truncated or corrupt compressed file ({exc})") from None


def _sniff_format(path: Path) -> str | None:
    with open_text(path) as fh:
        for line in fh:
            if line.strip():
                if line[0] == "@":
                    return "fastq"
                if line[0] == ">":
                    return "fasta"
                raise RasherError(f"{path} is neither fasta nor fastq (first line: {line.strip()[:40]!r})")
    return None


def read_records(path: Path) -> Iterator[Record]:
    """Iterate over the records of a fasta (multi-line allowed) or fastq (4 lines per record) file.

    A truncated or corrupt compressed file is a RasherError (decompression errors are not). Reading stops when
    RASHER is interrupted (a sample's reads can take minutes to go through)."""
    try:
        for count, rec in enumerate(_read_records(path)):
            if count % 4096 == 0 and stopping():
                raise RasherError(f"{path}: reading stopped: RASHER is stopping")
            yield rec
    except DECOMPRESSION_ERRORS as exc:
        raise RasherError(f"{path}: truncated or corrupt compressed file ({exc})") from None


def _read_records(path: Path) -> Iterator[Record]:
    fmt = sniff_format(path)
    if fmt is None:
        return
    with open_text(path) as fh:
        if fmt == "fasta":
            header, chunks = None, []
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None:
                        yield Record(header, "".join(chunks))
                    header, chunks = line[1:], []
                elif line and header is not None:
                    chunks.append(line)
            if header is not None:
                yield Record(header, "".join(chunks))
            return
        while True:
            header = fh.readline()
            if not header:
                return
            if not header.strip():
                continue
            seq, plus, qual = fh.readline().rstrip(), fh.readline(), fh.readline().rstrip()
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise RasherError(f"{path}: malformed or truncated fastq record {header.strip()[:60]!r}")
            yield Record(header[1:].rstrip(), seq, qual)


@dataclass
class ReadStats:
    """Number of reads, bases and (optionally) read lengths, for N50."""

    reads: int = 0
    bases: int = 0
    lengths: list[int] = field(default_factory=list, repr=False)

    def add(self, length: int, keep_length: bool = True) -> None:
        self.reads += 1
        self.bases += length
        if keep_length:
            self.lengths.append(length)

    @property
    def n50(self) -> int:
        return n50(self.lengths)


def n50(lengths: Iterable[int]) -> int:
    ordered = sorted(lengths, reverse=True)
    half, running = sum(ordered) / 2, 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return 0


def read_stats(paths: Iterable[Path]) -> ReadStats:
    stats = ReadStats()
    for path in paths:
        for rec in read_records(path):
            stats.add(len(rec.seq))
    return stats


def sequence_lengths(path: Path) -> list[tuple[str, int]]:
    return [(rec.name, len(rec.seq)) for rec in read_records(path)]


def check_reference(path: Path) -> list[tuple[str, int]]:
    """Validate a fasta reference and return its (name, length) pairs."""
    if not path.is_file():
        raise RasherError(f"Reference file not found: {path}")
    if sniff_format(path) != "fasta":
        raise RasherError(f"The reference must be a fasta file: {path}")
    lengths = sequence_lengths(path)
    if not lengths or sum(n for _, n in lengths) == 0:
        raise RasherError(f"The reference contains no sequence: {path}")
    names = [name for name, _ in lengths]
    if len(set(names)) != len(names):
        raise RasherError(f"The reference has duplicate sequence names: {path}")
    return lengths


_NOT_ACGT = re.compile(r"[^ACGT]")


def acgtn(seq: str) -> str:
    """Upper case, and N for anything that is not A, C, G or T (IUPAC codes, soft-masked bases)."""
    return _NOT_ACGT.sub("N", seq.upper())


def write_fasta(path: Path, records: Iterable[Record], width: int = 80) -> int:
    """Write records as plain fasta wrapped at `width`. Returns the number of records written."""
    count = 0
    with open(path, "w", encoding="ascii", errors="replace") as out:
        for rec in records:
            out.write(f">{rec.header}\n")
            for i in range(0, len(rec.seq), width):
                out.write(rec.seq[i: i + width] + "\n")
            count += 1
    return count


def canonical_name(name: str) -> str:
    """A read's name without a trailing /1 or /2: the mate suffixes that aligners strip (SAM spec), so a name
    reported by an aligner matches the read whatever the fastq calls it."""
    return name[:-2] if name.endswith(("/1", "/2")) else name


def extract_reads(inputs: Iterable[Path], names: set[str], output: Path) -> tuple[ReadStats, ReadStats]:
    """Copy the reads whose canonical name is in `names` (canonical names themselves) from `inputs` to `output`.

    Returns statistics for all input reads (without lengths) and for the extracted reads.
    """
    total, kept = ReadStats(), ReadStats()
    with open_write(output) as out:
        for path in inputs:
            for rec in read_records(path):
                total.add(len(rec.seq), keep_length=False)
                if canonical_name(rec.name) in names:
                    kept.add(len(rec.seq))
                    out.write(rec.format())
    return total, kept


def extract_pairs(pairs: Iterable[tuple[Path, Path]], names: set[str], out1: Path,
                  out2: Path) -> tuple[ReadStats, ReadStats]:
    """Copy the read pairs of which either mate's canonical name is in `names` (canonical names themselves) to
    `out1`/`out2`: a pair stays whole even when only one mate matched the reference.

    Each (R1, R2) file pair is read in lockstep: mates must come in the same order, as sequencers and the SRA
    write them (checked; identical mate names and `/1 /2` names both pass). Statistics count both mates.
    """
    total, kept = ReadStats(), ReadStats()
    with open_write(out1) as o1, open_write(out2) as o2:
        for r1_path, r2_path in pairs:
            for i, (rec1, rec2) in enumerate(zip_longest(read_records(r1_path), read_records(r2_path)), 1):
                if rec1 is None or rec2 is None:
                    raise RasherError(f"{r1_path} and {r2_path} do not hold the same number of reads: "
                                      "not mates of each other")
                name1, name2 = canonical_name(rec1.name), canonical_name(rec2.name)
                if name1 != name2:
                    raise RasherError(f"{r1_path} and {r2_path} desynchronise at read {i}: {rec1.name!r} vs "
                                      f"{rec2.name!r} (R1 and R2 must list the mates in the same order)")
                total.add(len(rec1.seq), keep_length=False)
                total.add(len(rec2.seq), keep_length=False)
                if name1 in names:
                    kept.add(len(rec1.seq))
                    kept.add(len(rec2.seq))
                    o1.write(rec1.format())
                    o2.write(rec2.format())
    return total, kept


def subsample_reads(in1: Path, in2: Path | None, out1: Path, out2: Path | None,
                    target_bases: int) -> ReadStats:
    """Keep an evenly spread subset of the reads, totalling about `target_bases`, pairs kept whole.

    Systematic sampling of the units (a pair, or a single read): unit i is kept when floor((i+1)*f) exceeds
    floor(i*f), f being the fraction of the bases to keep. Deterministic (no seed), order-preserving, and
    even across the file, so a tile- or time-ordered run stays represented. With `in2`, reads are paired and
    read in lockstep (not checked here: the baiting wrote these files pair by pair)."""
    bases = sum(len(rec.seq) for rec in read_records(in1))
    if in2 is not None:
        bases += sum(len(rec.seq) for rec in read_records(in2))
    fraction = min(1.0, target_bases / bases) if bases else 1.0
    kept = ReadStats()

    def keep(i: int) -> bool:
        return int((i + 1) * fraction) > int(i * fraction)

    second = read_records(in2) if in2 is not None else None
    with open_write(out1) as o1, (open_write(out2) if out2 is not None else nullcontext()) as o2:
        for i, rec1 in enumerate(read_records(in1)):
            rec2 = next(second, None) if second is not None else None
            if not keep(i):
                continue
            kept.add(len(rec1.seq))
            o1.write(rec1.format())
            if rec2 is not None:
                kept.add(len(rec2.seq))
                o2.write(rec2.format())
    return kept


def concatenate(inputs: Iterable[Path], output: Path) -> None:
    with open_write(output) as out:
        for path in inputs:
            for rec in read_records(path):
                out.write(rec.format())
