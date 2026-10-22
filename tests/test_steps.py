import shutil

import pytest

from rasher.steps import _fill_from, count_circular


def test_count_circular_headers(tmp_path):
    p = tmp_path / "a.fasta"
    p.write_text(">u1ctg_len-70000_circular-yes_depth-30 mult=1\nA\n"   # myloasm
                 ">u2ctg_len-9000_circular-no_depth-30\nA\n"
                 ">Utg30 LN:i:69995 RC:i:15 XO:i:1\nA\n"                  # Raven
                 ">Utg31 LN:i:100 XO:i:0\nA\n"
                 ">contig_1 circular=true\nA\n>plain\nA\n")
    assert count_circular(p) == 3


def test_fill_from_reference_only_at_the_ends():
    assert _fill_from("NNACGTNNTN", "GGACGTAAAA") == "GGACGTNNTA"
    assert _fill_from("ACGT", "") == "ACGT"
    assert _fill_from("NNNN", "ACGT") == "NNNN"


def test_publish_marks_circular_contigs(tmp_path):
    from rasher.seqio import read_records
    from rasher.steps import AssemblyDirs, _publish_assembly
    raw = tmp_path / "raw.fasta"
    raw.write_text(">contig_1\nACGT\n>u2ctg_circular-yes_depth-3\nAC\n>contig_3\nGG\n")
    dirs = AssemblyDirs(tmp_path / "asm")
    _publish_assembly("s1", raw, dirs, rename=False, circular={"contig_3"})
    headers = [r.header for r in read_records(dirs.assemblies / "s1.fasta")]
    assert headers == ["s1_contig_1", "s1_u2ctg_circular-yes_depth-3 circular=true", "s1_contig_3 circular=true"]


def test_pipe_failure_names_the_program_that_failed(tmp_path):
    import pytest

    from rasher.tools import ToolError, run
    with pytest.raises(ToolError, match="^false failed with exit code 1"):
        run([["yes"], ["false"]], tmp_path / "log")
    with pytest.raises(ToolError, match="^sh failed with exit code 7"):
        run([["sh", "-c", "exit 7"], ["cat"]], tmp_path / "log")


def test_published_assemblies_are_upper_case_acgtn(tmp_path):
    from rasher.seqio import read_records
    from rasher.steps import AssemblyDirs, _publish_assembly
    raw = tmp_path / "raw.fasta"
    raw.write_text(">c1\nacgtRYKMn\n")
    dirs = AssemblyDirs(tmp_path / "asm")
    _publish_assembly("s1", raw, dirs, rename=False)
    assert next(read_records(dirs.assemblies / "s1.fasta")).seq == "ACGTNNNNN"


def test_no_program_starts_once_rasher_is_stopping(tmp_path):
    import pytest

    from rasher.tools import ToolError, allow_programs, kill_running, run
    kill_running()
    try:
        with pytest.raises(ToolError, match="not started: RASHER is stopping"):
            run(["true"], tmp_path / "log")
    finally:
        allow_programs()
    run(["true"], tmp_path / "log")


def test_a_second_interruption_during_the_kill_leaves_no_program_running(tmp_path, monkeypatch):
    import signal
    import subprocess
    import threading

    import pytest

    from rasher import tools
    real_kill, seen = tools._kill, []

    def kill(proc):  # A second Ctrl-C arrives once the first program is killed, before the second one is
        seen.append(proc)
        if len(seen) == 2:
            raise KeyboardInterrupt
        real_kill(proc)

    monkeypatch.setattr(tools, "_kill", kill)
    timer = threading.Timer(0.5, signal.pthread_kill, (threading.main_thread().ident, signal.SIGINT))
    timer.start()  # The first Ctrl-C
    try:
        with pytest.raises(KeyboardInterrupt):
            tools.run([["sleep", "30"], ["sleep", "30"]], tmp_path / "log")
    finally:
        timer.cancel()
        monkeypatch.setattr(tools, "_kill", real_kill)
    assert len(seen) == 2
    tools.kill_running()  # As RASHER does at the end of an interrupted run
    for proc in seen:
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            real_kill(proc)
            raise AssertionError(f"{proc.args} still running") from None


def test_templated_assembly_notes_reference_sequences_without_reads(tmp_path, monkeypatch):
    from pathlib import Path

    import rasher.steps as steps

    def fake_run(cmd, log_file, **kwargs):  # The consensus has only the reference sequence with reads
        if "consensus" in cmd:
            Path(cmd[cmd.index("-o") + 1]).write_text(">chr\nACGTACGTAC\n")

    monkeypatch.setattr(steps, "run", fake_run)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr\nACGTACGTAC\n>plasmid\nGGGGCCCC\n")
    res = steps.assemble_samtools("s1", tmp_path / "reads.fq", None, reference, steps.AssemblyDirs(tmp_path / "asm"),
                                  tmp_path / "logs", fill_gaps=False, threads=1, reference_length=18)
    assert "no reads on reference sequence(s) plasmid: not in the assembly" in res.notes


def test_only_a_sample_folder_in_the_step_folder_is_deleted(tmp_path):
    import pytest

    from rasher.steps import AssemblyDirs, SampleFailed
    root = tmp_path / "out" / "3_assembled"
    (root / "s1").mkdir(parents=True)
    (root / "s1" / "old.fasta").write_text(">x\nA\n")
    (tmp_path / "out" / "reads.fastq").write_text("@r\nA\n+\nI\n")
    dirs = AssemblyDirs(root)
    assert dirs.clear_work("s1") == root / "s1" and not (root / "s1").exists()
    for name in ("..", ".", "s1/.."):
        with pytest.raises(SampleFailed, match="refusing to delete"):
            dirs.clear_work(name)
    assert (tmp_path / "out" / "reads.fastq").is_file() and root.is_dir()


def test_logged_command_line_can_be_run_again(tmp_path):
    """The working folder of a program is in its logged command line."""
    from rasher.tools import run
    work = tmp_path / "my work"
    work.mkdir()
    run(["sh", "-c", "pwd"], tmp_path / "log", cwd=work)
    run(["echo", "a b"], tmp_path / "log", stdout=tmp_path / "out put.txt")
    lines = (tmp_path / "log").read_text().splitlines()
    assert lines[:2] == [f"$ (cd '{work}' && sh -c pwd)", str(work)]
    assert lines[2] == f"$ echo 'a b' > '{tmp_path / 'out put.txt'}'"


@pytest.mark.skipif(not shutil.which("bbduk.sh"), reason="BBDuk is not installed")
@pytest.mark.parametrize("names", [["v.FASTQ.GZ"], ["y.FQ.GZ"], ["z.fastq"], ["x.fastq.gz"], ["w.fq.gz"],
                                   ["a.FQ.GZ", "b.fastq"]])
def test_bbduk_reads_files_whatever_their_extension_says(tmp_path, names):
    """BBTools picks the format and compression from the extension (case-sensitively); RASHER from the content."""
    import gzip
    import random

    from rasher.samples import Sample
    from rasher.steps import bait_bbduk
    rng = random.Random(1)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">r\n{genome}\n")
    files = []
    for i, name in enumerate(names):
        reads = "".join(f"@r{i}_{j}\n{genome[j * 100:j * 100 + 300]}\n+\n{'I' * 300}\n" for j in range(5))
        path = tmp_path / name
        # Gzipped or not, against what the name says, for the two lower-case cases
        compressed = name.upper().endswith(".GZ") != (name in ("z.fastq", "x.fastq.gz"))
        path.write_bytes(gzip.compress(reads.encode()) if compressed else reads.encode())
        files.append(path)
    result = bait_bbduk(Sample("s", files, "fastq"), ref, tmp_path / "out", tmp_path, threads=1, kmer=31,
                        hdist=0, memory_gb=1)
    assert result.stats["Baited_reads"] == 5 * len(names)
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["s.fastq.gz"]  # No link left behind


@pytest.mark.skipif(not shutil.which("bbduk.sh"), reason="BBDuk is not installed")
def test_bbduk_keeps_a_pair_with_one_mate_matching(tmp_path):
    import random

    from rasher.samples import Sample
    from rasher.seqio import read_records
    from rasher.steps import bait_bbduk
    rng = random.Random(2)
    genome, other = ("".join(rng.choice("ACGT") for _ in range(3000)) for _ in range(2))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">r\n{genome}\n")
    pairs = []
    for lane, (first, second) in (("a", (genome, other)), ("b", (other, other))):  # Lane b: nothing matches
        r1 = _write_reads(tmp_path / f"{lane}_R1.fq", [(f"{lane}{j}/1", other[j * 100:j * 100 + 150])
                                                       for j in range(5)], "fastq")
        r2 = _write_reads(tmp_path / f"{lane}_R2.fq", [(f"{lane}{j}/2", (first if j == 2 else second)
                                                        [j * 100:j * 100 + 150]) for j in range(5)], "fastq")
        pairs.append((r1, r2))
    sample = Sample("s", [f for pair in pairs for f in pair], "fastq", pairs)
    res = bait_bbduk(sample, ref, tmp_path / "out", tmp_path, threads=1, kmer=31, hdist=0, memory_gb=1)
    assert [r.name for r in read_records(res.output)] == ["a2/1"]  # Its mate matched
    assert [r.name for r in read_records(res.output2)] == ["a2/2"]
    assert res.stats["Raw_reads"] == 20 and res.stats["Baited_reads"] == 2
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["s_R1.fastq.gz", "s_R2.fastq.gz"]


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 is not installed")
def test_off_reference_share(tmp_path):
    import random

    from rasher.steps import off_reference
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(20000))
    other = "".join(rng.choice("ACGT") for _ in range(5000))
    (tmp_path / "ref.fa").write_text(f">ref\n{ref}\n")
    (tmp_path / "asm.fa").write_text(f">c1 circular=true\n{ref[:15000]}\n>c2\n{other}\n"
                                     f">c3\n{ref[15000:]}{other[:2000]}\n")
    pct = off_reference(tmp_path / "asm.fa", tmp_path / "ref.fa", tmp_path / "log", 1)
    assert 25 < pct < 27  # 7000 of 27000 bases
    (tmp_path / "same.fa").write_text(f">c1\n{ref}\n")
    assert off_reference(tmp_path / "same.fa", tmp_path / "ref.fa", tmp_path / "log", 1) < 1


def _write_reads(path, reads, fmt="fasta"):
    """`reads`: (name, sequence) pairs."""
    path.write_text("".join(f">{n}\n{s}\n" if fmt == "fasta" else f"@{n}\n{s}\n+\n{'I' * len(s)}\n"
                            for n, s in reads))
    return path


def test_bait_result_counts_both_mates_and_names_both_outputs(tmp_path):
    from rasher.seqio import ReadStats
    from rasher.steps import SampleFailed, _bait_result
    r1, r2 = tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz"
    res = _bait_result([r1, r2], ReadStats(10, 1_000), ReadStats(4, 400))
    assert (res.output, res.output2) == (r1, r2)
    assert res.stats == {"Raw_reads": 10, "Raw_bases": 1_000, "Baited_reads": 4, "Baited_bases": 400,
                         "Baited_pct": "40.000"}
    assert _bait_result([r1], None, ReadStats(1, 10)).output2 is None
    with pytest.raises(SampleFailed, match="no reads matched") as exc:
        _bait_result([r1, r2], ReadStats(10, 1_000), ReadStats())
    assert exc.value.stats["Raw_reads"] == 10  # What is known of the failed sample


def _fake_minimap2(hits, calls):
    """A `run` reporting the reads of `hits` (names as the aligner gives them) as aligned, in each file."""
    from pathlib import Path

    def run(cmds, log_file, *, stdout=None, **kwargs):
        calls.append(cmds)
        if stdout is not None:
            names = {line[1:].split()[0] for line in Path(cmds[-1]).read_text().splitlines()[::4]}
            Path(stdout).write_text("".join(f"{n}\t150\t0\t150\t+\tref\n" for n in hits if n in names))
    return run


def test_minimap2_bait_keeps_a_pair_with_one_mate_aligned(tmp_path, monkeypatch):
    import rasher.steps as steps
    from rasher.samples import Sample
    from rasher.seqio import read_records
    pairs = []
    for lane in ("a", "b"):
        r1 = _write_reads(tmp_path / f"{lane}_R1.fq", [(f"{lane}{i}/1", "ACGT" * 10) for i in range(3)], "fastq")
        r2 = _write_reads(tmp_path / f"{lane}_R2.fq", [(f"{lane}{i}/2", "TTGCA" * 10) for i in range(3)], "fastq")
        pairs.append((r1, r2))
    calls = []
    monkeypatch.setattr(steps, "run", _fake_minimap2({"a1/1", "b2/2"}, calls))  # One mate of each pair
    sample = Sample("s", [f for pair in pairs for f in pair], "fastq", pairs)
    res = steps.bait_minimap2(sample, tmp_path / "ref.fa", tmp_path / "out", tmp_path / "logs", 4, keep_bam=True)
    out = tmp_path / "out"
    assert (res.output, res.output2) == (out / "s_R1.fastq.gz", out / "s_R2.fastq.gz")
    assert [r.name for r in read_records(res.output)] == ["a1/1", "b2/1"]
    assert [r.name for r in read_records(res.output2)] == ["a1/2", "b2/2"]
    assert res.stats["Raw_reads"] == 12 and res.stats["Baited_reads"] == 4
    # Each file on its own (minimap2 -x sr would pair two files, and refuses more); the BAM from the baited pairs
    assert [c[:3] + c[-1:] for c in calls[:4]] == [["minimap2", "-x", "sr", str(f)] for f in sample.files]
    assert calls[4][0][:3] == ["minimap2", "-ax", "sr"] and calls[4][0][-2:] == [str(res.output), str(res.output2)]
    assert calls[5] == ["samtools", "index", str(out / "s.bam")]
    assert sorted(p.name for p in out.iterdir()) == ["s_R1.fastq.gz", "s_R2.fastq.gz"]  # No PAF left behind


def test_minimap2_bait_of_single_reads_matches_the_names_without_mate_suffix(tmp_path, monkeypatch):
    import rasher.steps as steps
    from rasher.samples import Sample
    from rasher.seqio import read_records
    reads = _write_reads(tmp_path / "r.fq", [("x/1", "ACGT" * 10), ("y/1", "ACGT" * 10)], "fastq")
    monkeypatch.setattr(steps, "run", _fake_minimap2({"x/1"}, []))
    res = steps.bait_minimap2(Sample("s", [reads], "fastq"), tmp_path / "ref.fa", tmp_path / "out", tmp_path, 1,
                              keep_bam=False)
    assert res.output == tmp_path / "out" / "s.fastq.gz" and res.output2 is None
    assert [r.name for r in read_records(res.output)] == ["x/1"]
    monkeypatch.setattr(steps, "run", _fake_minimap2(set(), []))
    with pytest.raises(steps.SampleFailed, match="no reads matched"):
        steps.bait_minimap2(Sample("s", [reads], "fastq"), tmp_path / "ref.fa", tmp_path / "out", tmp_path, 1,
                            keep_bam=False)
    assert not res.output.exists()  # A failed sample leaves no reads


def test_fasta_pairs_are_filtered_whole_and_subsampled(tmp_path):
    from rasher.seqio import read_records
    from rasher.steps import filter_fastp
    r1 = _write_reads(tmp_path / "in_R1.fasta", [(f"p{i}/1", "A" * 100) for i in range(10)])
    r2 = _write_reads(tmp_path / "in_R2.fasta", [(f"p{i}/2", "C" * (40 if i == 3 else 100)) for i in range(10)])
    out = tmp_path / "out"
    res = filter_fastp("s", r1, r2, "fasta", out, tmp_path / "logs", min_length=50, target_bases=10_000,
                       threads=1)
    assert (res.output, res.output2) == (out / "s_R1.fasta.gz", out / "s_R2.fasta.gz")
    assert [r.name for r in read_records(res.output)] == [f"p{i}/1" for i in range(10) if i != 3]  # p3: short R2
    assert [r.name for r in read_records(res.output2)] == [f"p{i}/2" for i in range(10) if i != 3]
    assert res.stats == {"Filtered_reads": 18, "Filtered_bases": 1_800, "Filtered_N50": 100}
    res = filter_fastp("s", r1, r2, "fasta", out, tmp_path / "logs", min_length=50, target_bases=900, threads=1)
    assert res.stats["Filtered_reads"] == 8  # Half the bases: every other pair, both mates
    assert sorted(p.name for p in out.iterdir()) == ["s_R1.fasta.gz", "s_R2.fasta.gz"]  # No temporary file left


def test_filtering_single_fasta_reads(tmp_path):
    from rasher.seqio import read_records
    from rasher.steps import filter_fastp
    reads = _write_reads(tmp_path / "in.fasta", [("a", "A" * 100), ("b", "A" * 30), ("c", "A" * 60)])
    res = filter_fastp("s", reads, None, "fasta", tmp_path / "out", tmp_path / "logs", min_length=50,
                       target_bases=10_000, threads=1)
    assert res.output == tmp_path / "out" / "s.fasta.gz" and res.output2 is None
    assert [r.name for r in read_records(res.output)] == ["a", "c"]


def test_filtering_fails_when_no_read_is_left(tmp_path):
    from rasher.steps import SampleFailed, filter_fastp
    reads = _write_reads(tmp_path / "in.fasta", [("a", "ACGT"), ("b", "ACGTACGT")])
    with pytest.raises(SampleFailed, match=r"no reads left after filtering \(all shorter than 100 bp\?\)"):
        filter_fastp("s", reads, None, "fasta", tmp_path / "out", tmp_path / "logs", min_length=100,
                     target_bases=1_000, threads=1)
    assert list((tmp_path / "out").iterdir()) == []  # Neither an empty output nor a temporary file


def _fake_fastp(calls):
    """A `run` for which fastp keeps every read."""
    import shutil
    from pathlib import Path

    def run(cmd, log_file, **kwargs):
        calls.append(cmd)
        for given, written in (("--in1", "--out1"), ("--in2", "--out2")):
            if given in cmd:
                shutil.copyfile(cmd[cmd.index(given) + 1], cmd[cmd.index(written) + 1])
        Path(cmd[cmd.index("--json") + 1]).write_text("{}")
    return run


def test_fastp_filters_the_pairs_into_hidden_files_and_reports_to_the_log_folder(tmp_path, monkeypatch):
    import rasher.steps as steps
    calls = []
    monkeypatch.setattr(steps, "run", _fake_fastp(calls))
    r1 = _write_reads(tmp_path / "in_R1.fq", [("p/1", "A" * 100)], "fastq")
    r2 = _write_reads(tmp_path / "in_R2.fq", [("p/2", "C" * 100)], "fastq")
    out, logs = tmp_path / "out", tmp_path / "logs"
    logs.mkdir()
    res = steps.filter_fastp("s", r1, r2, "fastq", out, logs, min_length=50, target_bases=1_000, threads=3)
    assert calls == [["fastp", "--in1", str(r1), "--out1", str(out / ".s.trimmed_R1.fastq.gz"),
                      "--in2", str(r2), "--out2", str(out / ".s.trimmed_R2.fastq.gz"), "--length_required", "50",
                      "--thread", "3", "--json", str(logs / "s.fastp.json"), "--html", str(logs / "s.fastp.html")]]
    assert (res.output, res.output2) == (out / "s_R1.fastq.gz", out / "s_R2.fastq.gz")
    assert res.stats["Filtered_reads"] == 2
    assert sorted(p.name for p in out.iterdir()) == ["s_R1.fastq.gz", "s_R2.fastq.gz"]
    calls.clear()
    res = steps.filter_fastp("t", r1, None, "fastq", out, logs, min_length=50, target_bases=1_000, threads=1)
    assert "--in2" not in calls[0] and calls[0][4] == str(out / ".t.trimmed.fastq.gz")
    assert res.output == out / "t.fastq.gz" and res.output2 is None
    assert not (out / ".t.trimmed.fastq.gz").exists()


def _fake_assembler(calls, contigs=">c1\nACGTACGT\n>c2\nGGCC\n"):
    """A `run` for which SPAdes and SKESA write `contigs` (and SPAdes a graph)."""
    from pathlib import Path

    def run(cmd, log_file, *, cwd=None, **kwargs):
        calls.append((cmd, cwd))
        if cmd[0] == "spades.py":
            work = Path(cmd[cmd.index("-o") + 1])
            work.mkdir(exist_ok=True)
            (work / "contigs.fasta").write_text(contigs)
            (work / "assembly_graph_with_scaffolds.gfa").write_text("S\t1\tACGT\n")
        elif cmd[0] == "skesa":
            assert all((cwd / name).resolve().is_file() for name in cmd[cmd.index("--reads") + 1].split(","))
            (cwd / cmd[cmd.index("--contigs_out") + 1]).write_text(contigs)
    return run


def test_spades_assembles_pairs_or_single_reads(tmp_path, monkeypatch):
    import rasher.steps as steps
    from rasher.seqio import read_records
    calls = []
    monkeypatch.setattr(steps, "run", _fake_assembler(calls))
    monkeypatch.setattr(steps, "which", lambda name: None)  # No Bandage picture
    r1, r2 = tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz"
    dirs = steps.AssemblyDirs(tmp_path / "asm")
    res = steps.assemble_spades("s", r1, r2, dirs, tmp_path / "logs", threads=2, memory_gb=8, reference_length=24)
    assert calls[0][0] == ["spades.py", "--isolate", "-1", str(r1), "-2", str(r2), "-t", "2", "-m", "8", "-o",
                           str(dirs.work("s"))]
    assert [r.header for r in read_records(res.output)] == ["s_c1", "s_c2"]
    assert (dirs.graphs / "s.gfa").is_file()
    assert res.stats == {"Contigs": 2, "Assembly_length": 12, "Largest_contig": 8, "Assembly_N50": 8,
                         "Length_vs_reference": "0.500", "Circular_contigs": 0, "Assembly_depth": "NA"}
    steps.assemble_spades("s", r1, None, dirs, tmp_path / "logs", threads=2, memory_gb=8, reference_length=24)
    assert calls[1][0][2:4] == ["-s", str(r1)] and "-2" not in calls[1][0]


def test_spades_works_through_links_when_a_path_has_a_space(tmp_path, monkeypatch):
    import rasher.steps as steps
    calls = []
    monkeypatch.setattr(steps, "run", _fake_assembler(calls))
    monkeypatch.setattr(steps, "which", lambda name: None)
    folder = tmp_path / "my run"
    folder.mkdir()
    r1 = _write_reads(folder / "s_R1.fastq", [("p/1", "ACGT")], "fastq")
    r2 = _write_reads(folder / "s_R2.fastq", [("p/2", "ACGT")], "fastq")
    dirs = steps.AssemblyDirs(folder / "asm")
    res = steps.assemble_spades("s", r1, r2, dirs, tmp_path / "logs", threads=1, memory_gb=1, reference_length=12)
    cmd = calls[0][0]
    assert not any(" " in arg for arg in cmd)
    assert [cmd[3].rsplit("/", 1)[1], cmd[5].rsplit("/", 1)[1]] == ["reads-1.fastq", "reads-2.fastq"]
    assert (dirs.work("s") / "contigs.fasta").is_file() and res.stats["Contigs"] == 2
    assert not steps.Path(cmd[-1]).parent.exists()  # The links' folder is gone


def test_spades_without_contigs_points_at_its_log(tmp_path, monkeypatch):
    import rasher.steps as steps
    monkeypatch.setattr(steps, "run", lambda *args, **kwargs: None)
    dirs = steps.AssemblyDirs(tmp_path / "asm")
    with pytest.raises(steps.SampleFailed, match="SPAdes produced no assembly .see .*spades.log"):
        steps.assemble_spades("s", tmp_path / "r.fq", None, dirs, tmp_path, threads=1, memory_gb=1,
                              reference_length=10)


def test_skesa_reads_pairs_through_comma_free_links(tmp_path, monkeypatch):
    import rasher.steps as steps
    calls = []
    monkeypatch.setattr(steps, "run", _fake_assembler(calls))
    folder = tmp_path / "a,b"  # A comma would split a path of SKESA's --reads
    folder.mkdir()
    r1 = _write_reads(folder / "s_R1.fastq", [("p/1", "ACGT")], "fastq")
    r2 = _write_reads(folder / "s_R2.fastq", [("p/2", "ACGT")], "fastq")
    dirs = steps.AssemblyDirs(tmp_path / "asm")
    dirs.graphs.mkdir(parents=True)
    (dirs.graphs / "s.gfa").write_text("S\t1\tA\n")  # From an earlier SPAdes assembly
    res = steps.assemble_skesa("s", r1, r2, dirs, tmp_path / "logs", threads=2, memory_gb=6, reference_length=24)
    work = dirs.work("s")
    assert calls == [(["skesa", "--reads", "reads_R1.fastq,reads_R2.fastq", "--contigs_out", "contigs.fasta",
                       "--cores", "2", "--memory", "6"], work)]
    assert sorted(p.name for p in work.iterdir()) == ["contigs.fasta"]  # The links are gone
    assert not (dirs.graphs / "s.gfa").exists()  # SKESA draws no graph: an older one would mislead
    assert res.stats["Contigs"] == 2 and res.stats["Assembly_depth"] == "NA"
    calls.clear()
    steps.assemble_skesa("s", r1, None, dirs, tmp_path / "logs", threads=2, memory_gb=6, reference_length=24)
    assert calls[0][0][1:3] == ["--reads", "reads.fastq"]


def test_skesa_without_contigs_points_at_the_log(tmp_path, monkeypatch):
    import rasher.steps as steps
    monkeypatch.setattr(steps, "run", _fake_assembler([], contigs=""))
    reads = _write_reads(tmp_path / "r.fq", [("p", "ACGT")], "fastq")
    with pytest.raises(steps.SampleFailed, match=r"SKESA produced no contigs \(see .*s\.log\)"):
        steps.assemble_skesa("s", reads, None, steps.AssemblyDirs(tmp_path / "asm"), tmp_path, threads=1,
                             memory_gb=1, reference_length=10)


def test_templated_assembly_aligns_pairs_as_pairs_with_the_illumina_profile(tmp_path, monkeypatch):
    from pathlib import Path

    import rasher.steps as steps
    calls = []

    def fake_run(cmd, log_file, **kwargs):
        calls.append(cmd)
        if "consensus" in cmd:
            Path(cmd[cmd.index("-o") + 1]).write_text(">chr\nACGT\n")

    monkeypatch.setattr(steps, "run", fake_run)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr\nACGT\n")
    r1, r2 = tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz"
    steps.assemble_samtools("s", r1, r2, reference, steps.AssemblyDirs(tmp_path / "asm"), tmp_path,
                            threads=1, reference_length=4, fill_gaps=False)
    assert calls[0][0][:3] == ["minimap2", "-ax", "sr"] and calls[0][0][-3:] == [str(reference), str(r1), str(r2)]
    assert calls[2][calls[2].index("-X") + 1] == "hiseq"
