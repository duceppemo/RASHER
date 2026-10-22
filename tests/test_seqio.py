import pytest

from rasher import RasherError
from rasher.seqio import (
    check_reference,
    extract_reads,
    n50,
    read_records,
    read_stats,
    sniff_format,
    split_extension,
)


@pytest.mark.parametrize("name, expected", [
    ("s1.fastq.gz", ("s1", "fastq")), ("s1.fq", ("s1", "fastq")), ("a.b.fasta", ("a.b", "fasta")),
    ("x.FNA.GZ", ("x", "fasta")), ("notes.txt", None), (".fastq", None), ("s.fastq.bz2", None),
])
def test_split_extension(name, expected):
    assert split_extension(name) == expected


def test_read_fastq_gz_and_plain(fastq):
    for gz in (True, False):
        path = fastq(f"r{gz}.fastq" + (".gz" if gz else ""), [("r1", "ACGT"), ("r2", "AAAAAA")], gz)
        recs = list(read_records(path))
        assert [(r.name, r.seq) for r in recs] == [("r1", "ACGT"), ("r2", "AAAAAA")]
        assert recs[0].header == "r1 extra"


def test_gzip_detected_by_content_not_name(tmp_path, fastq):
    path = fastq("mislabelled.fastq", [("r1", "ACGT")], gz=True)  # gzipped, no .gz
    assert [r.seq for r in read_records(path)] == ["ACGT"]


def test_multiline_fasta(tmp_path):
    p = tmp_path / "m.fasta"
    p.write_text(">a desc\nACG\nTT\n\n>b\nGG\n")
    assert [(r.name, r.seq) for r in read_records(p)] == [("a", "ACGTT"), ("b", "GG")]


def test_truncated_fastq_is_an_error(tmp_path):
    p = tmp_path / "t.fastq"
    p.write_text("@r1\nACGT\n+\nII\n")
    with pytest.raises(RasherError, match="malformed or truncated"):
        list(read_records(p))


def test_sniff_format(tmp_path):
    (tmp_path / "e.fastq").write_text("\n\n")
    assert sniff_format(tmp_path / "e.fastq") is None
    (tmp_path / "x.fasta").write_text("hello\n")
    with pytest.raises(RasherError, match="neither fasta nor fastq"):
        sniff_format(tmp_path / "x.fasta")


def test_n50():
    assert n50([]) == 0
    assert n50([10]) == 10
    assert n50([2, 3, 4, 10]) == 10  # 10 >= 19 / 2
    assert n50([5, 5, 5, 5]) == 5
    assert n50([1, 2, 3, 4, 5]) == 4  # 5 + 4 = 9 >= 7.5


def test_read_stats(fastq):
    stats = read_stats([fastq("a.fastq.gz", [("r1", "ACGT"), ("r2", "A" * 10)])])
    assert (stats.reads, stats.bases, stats.n50) == (2, 14, 10)


def test_extract_reads(fastq, tmp_path):
    src1 = fastq("a.fastq.gz", [("r1", "ACGT"), ("r2", "GGGG")])
    src2 = fastq("b.fastq", [("r3", "TTTTTT")], gz=False)
    out = tmp_path / "out.fastq.gz"
    total, kept = extract_reads([src1, src2], {"r2", "r3", "absent"}, out)
    assert (total.reads, total.bases, kept.reads, kept.bases) == (3, 14, 2, 10)
    assert [r.name for r in read_records(out)] == ["r2", "r3"]
    assert next(read_records(out)).qual == "IIII"


def test_check_reference(fasta):
    assert check_reference(fasta("ok.fasta", [("a", "ACGT"), ("b", "GG")])) == [("a", 4), ("b", 2)]
    with pytest.raises(RasherError, match="duplicate"):
        check_reference(fasta("dup.fasta", [("a", "ACGT"), ("a", "GG")]))
    with pytest.raises(RasherError, match="no sequence"):
        check_reference(fasta("empty.fasta", [("a", "")]))


def test_check_reference_rejects_fastq_and_missing(fastq, tmp_path):
    with pytest.raises(RasherError, match="must be a fasta"):
        check_reference(fastq("r.fastq.gz", [("r1", "ACGT")]))
    with pytest.raises(RasherError, match="not found"):
        check_reference(tmp_path / "missing.fasta")


def test_truncated_or_corrupt_gzip_is_a_rasher_error(tmp_path, fastq):
    import random
    rng = random.Random(3)  # Random reads: a compressed file large enough to corrupt its middle
    good = fastq("g.fastq.gz", [(f"r{i}", "".join(rng.choice("ACGT") for _ in range(800))) for i in range(300)])
    data = good.read_bytes()
    (tmp_path / "t.fastq.gz").write_bytes(data[: len(data) // 2])
    corrupt = bytearray(data)
    corrupt[10:60] = bytes(50)  # Right after the gzip header: an invalid deflate block
    (tmp_path / "c.fastq.gz").write_bytes(bytes(corrupt))
    for name in ("t.fastq.gz", "c.fastq.gz"):
        with pytest.raises(RasherError, match="truncated or corrupt compressed file"):
            list(read_records(tmp_path / name))


def test_non_ascii_headers_do_not_break_writing(tmp_path):
    from rasher.seqio import Record, open_write, write_fasta
    write_fasta(tmp_path / "a.fasta", [Record("Solanum × sp", "ACGT")])
    with open_write(tmp_path / "b.fastq.gz") as fh:
        fh.write(Record("é", "AC", "II").format())
    assert (tmp_path / "a.fasta").read_text().startswith(">Solanum ? sp")


def test_acgtn():
    from rasher.seqio import acgtn
    assert acgtn("acgtRYN-*x") == "ACGTNNNNNN"


def test_reading_stops_when_rasher_is_interrupted(fastq):
    import pytest

    from rasher import RasherError, tools
    from rasher.seqio import read_records
    path = fastq("r.fq.gz", [(f"r{i}", "ACGT") for i in range(10)])
    tools.kill_running()
    with pytest.raises(RasherError, match="RASHER is stopping"):
        list(read_records(path))
    tools.allow_programs()
    assert len(list(read_records(path))) == 10


def test_a_utf8_bom_is_skipped(tmp_path):
    import gzip
    bom = b"\xef\xbb\xbf"
    (tmp_path / "r.fasta").write_bytes(bom + b">a\nACGT\n")
    (tmp_path / "r.fastq.gz").write_bytes(gzip.compress(bom + b"@r1\nACGT\n+\nIIII\n"))
    assert sniff_format(tmp_path / "r.fasta") == "fasta" and sniff_format(tmp_path / "r.fastq.gz") == "fastq"
    assert [(r.name, r.seq) for r in read_records(tmp_path / "r.fasta")] == [("a", "ACGT")]
    assert [(r.name, r.qual) for r in read_records(tmp_path / "r.fastq.gz")] == [("r1", "IIII")]
    assert check_reference(tmp_path / "r.fasta") == [("a", 4)]
    (tmp_path / "x.fasta").write_bytes(b"\xef\xbb>a\nACGT\n")  # Not a BOM: still not fasta
    with pytest.raises(RasherError, match="neither fasta nor fastq"):
        sniff_format(tmp_path / "x.fasta")


# ----- Paired extraction and subsampling ------------------------------------------------------------------------

def test_canonical_name():
    from rasher.seqio import canonical_name
    assert canonical_name("read/1") == "read" and canonical_name("read/2") == "read"
    assert canonical_name("read") == "read" and canonical_name("ERR1.5") == "ERR1.5"


def test_extract_pairs_keeps_a_pair_when_either_mate_hits(fastq, tmp_path):
    from rasher.seqio import extract_pairs
    r1 = fastq("r1.fastq.gz", [("a", "ACGT"), ("b", "CCCC"), ("c", "GGGG")])
    r2 = fastq("r2.fastq.gz", [("a", "TTTT"), ("b", "AAAA"), ("c", "ACGT")])
    out1, out2 = tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz"
    total, kept = extract_pairs([(r1, r2)], {"a", "c"}, out1, out2)
    assert (total.reads, kept.reads) == (6, 4)
    from rasher.seqio import read_records
    assert [r.name for r in read_records(out1)] == ["a", "c"]
    assert [r.name for r in read_records(out2)] == ["a", "c"]


def test_extract_pairs_with_slash_mate_names(fastq, tmp_path):
    from rasher.seqio import extract_pairs, read_records
    r1 = fastq("r1.fastq.gz", [("a/1", "ACGT"), ("b/1", "CCCC")])
    r2 = fastq("r2.fastq.gz", [("a/2", "TTTT"), ("b/2", "AAAA")])
    out1, out2 = tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz"
    total, kept = extract_pairs([(r1, r2)], {"b"}, out1, out2)  # Aligners report "b", not "b/1"
    assert (total.reads, kept.reads) == (4, 2)
    assert [r.name for r in read_records(out1)] == ["b/1"]  # The fastq names are kept as they were
    assert [r.name for r in read_records(out2)] == ["b/2"]


def test_extract_pairs_desync_is_an_error(fastq, tmp_path):
    from rasher.seqio import extract_pairs
    r1 = fastq("r1.fastq.gz", [("a", "ACGT"), ("b", "CCCC")])
    r2 = fastq("r2.fastq.gz", [("b", "TTTT"), ("a", "AAAA")])
    with pytest.raises(RasherError, match="desynchronise at read 1"):
        extract_pairs([(r1, r2)], {"a"}, tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz")


@pytest.mark.parametrize("name1, name2, expected", [
    ("a", "a", True), ("a/1", "a/2", True), ("a/1", "a", True), ("SRR1.0.1", "SRR1.0.2", True),
    ("SRR1.0.2", "SRR1.0.1", True), ("x_1", "x_2", True), ("a", "b", False), ("SRR1.0.1", "SRR1.1.2", False),
    ("SRR1.1", "SRR1.3", False), ("SRR1.0.1", "SRR1.0.1/2", True), ("a.12", "a.21", False), ("x1", "y2", False),
    ("a.1", "a.3", False),
])
def test_mates_match(name1, name2, expected):
    from rasher.seqio import mates_match
    assert mates_match(name1, name2) is expected


def test_extract_pairs_with_fastq_dump_readids_names(fastq, tmp_path):
    """fastq-dump --readids: SRR1.<spot>.1 and SRR1.<spot>.2, which aligners report as they are."""
    from rasher.seqio import extract_pairs, read_records
    r1 = fastq("r1.fastq.gz", [(f"SRR1.{i}.1", "ACGT") for i in range(3)])
    r2 = fastq("r2.fastq.gz", [(f"SRR1.{i}.2", "TTTT") for i in range(3)])
    out1, out2 = tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz"
    total, kept = extract_pairs([(r1, r2)], {"SRR1.0.1", "SRR1.2.2"}, out1, out2)  # A hit on either mate
    assert (total.reads, kept.reads) == (6, 4)
    assert [r.name for r in read_records(out1)] == ["SRR1.0.1", "SRR1.2.1"]
    assert [r.name for r in read_records(out2)] == ["SRR1.0.2", "SRR1.2.2"]


def test_verify_pairs(fastq):
    from rasher.seqio import verify_pairs
    r1 = fastq("r1.fastq.gz", [("a/1", "ACGT"), ("b/1", "CC")])
    stats = verify_pairs(r1, fastq("r2.fastq.gz", [("a/2", "TTTTTT"), ("b/2", "A")]))
    assert (stats.reads, stats.bases) == (4, 13)
    with pytest.raises(RasherError, match="desynchronise at read 1: 'a/1' vs 'b/2'"):
        verify_pairs(r1, fastq("shuffled.fastq.gz", [("b/2", "A"), ("a/2", "TTTTTT")]))
    with pytest.raises(RasherError, match="same number of reads"):
        verify_pairs(r1, fastq("short.fastq.gz", [("a/2", "TTTTTT")]))


def test_extract_pairs_unequal_counts_is_an_error(fastq, tmp_path):
    from rasher.seqio import extract_pairs
    r1 = fastq("r1.fastq.gz", [("a", "ACGT"), ("b", "CCCC")])
    r2 = fastq("r2.fastq.gz", [("a", "TTTT")])
    with pytest.raises(RasherError, match="same number of reads"):
        extract_pairs([(r1, r2)], {"a"}, tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz")


def test_extract_reads_canonicalises_names(fastq, tmp_path):
    from rasher.seqio import extract_reads, read_records
    reads = fastq("r.fastq.gz", [("a/1", "ACGT"), ("b/1", "CCCC")])
    out = tmp_path / "o.fastq.gz"
    total, kept = extract_reads([reads], {"a"}, out)
    assert (total.reads, kept.reads) == (2, 1) and [r.name for r in read_records(out)] == ["a/1"]
    # Both mates of an interleaved file share their canonical name: kept together when either aligned
    interleaved = fastq("i.fastq.gz", [("a/1", "ACGT"), ("a/2", "TTTT"), ("b/1", "CCCC"), ("b/2", "GGGG")])
    extract_reads([interleaved], {"b"}, out)
    assert [r.name for r in read_records(out)] == ["b/1", "b/2"]


def test_subsample_pairs_whole_and_deterministic(fastq, tmp_path):
    from rasher.seqio import read_records, subsample_reads
    pairs = [(f"r{i}", "A" * 100) for i in range(100)]
    r1, r2 = fastq("r1.fastq.gz", pairs), fastq("r2.fastq.gz", pairs)
    out1, out2 = tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz"
    kept = subsample_reads(r1, r2, out1, out2, target_bases=5000)  # A quarter of 20,000 bases
    names1 = [r.name for r in read_records(out1)]
    assert kept.reads == 2 * len(names1) and 24 <= len(names1) <= 26  # Pairs whole, about a quarter
    assert names1 == [r.name for r in read_records(out2)]  # The same pairs in both outputs
    subsample_reads(r1, r2, tmp_path / "p1.fastq.gz", tmp_path / "p2.fastq.gz", target_bases=5000)
    assert (tmp_path / "p1.fastq.gz").read_bytes() != b"" and \
        [r.name for r in read_records(tmp_path / "p1.fastq.gz")] == names1  # Deterministic


def test_subsample_keeps_everything_when_under_target(fastq, tmp_path):
    from rasher.seqio import read_records, subsample_reads
    r1 = fastq("r1.fastq.gz", [("a", "ACGT"), ("b", "CCCC")])
    out = tmp_path / "o.fastq.gz"
    kept = subsample_reads(r1, None, out, None, target_bases=1_000_000)
    assert kept.reads == 2 and [r.name for r in read_records(out)] == ["a", "b"]


def test_subsample_single_end_fasta(fasta, tmp_path):
    from rasher.seqio import read_records, subsample_reads
    r = fasta("r.fasta", [(f"r{i}", "ACGT" * 25) for i in range(10)])
    out = tmp_path / "o.fasta"
    kept = subsample_reads(r, None, out, None, target_bases=300)
    assert kept.reads == 3 and len(list(read_records(out))) == 3
