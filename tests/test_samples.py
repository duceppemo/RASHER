import pytest

from rasher import RasherError
from rasher.samples import discover, read_sample_sheet


def test_single_file(fastq):
    samples = discover(fastq("s1.fastq.gz", [("r", "ACGT")]))
    assert [(s.name, s.fmt) for s in samples] == [("s1", "fastq")]


def test_folder_files_and_subfolders(tmp_path, fastq, fasta):
    d = tmp_path / "in"
    (d / "barcode01").mkdir(parents=True)
    (d / "unclassified").mkdir()
    (d / "notes").mkdir()
    (d / "notes" / "readme.txt").write_text("hi")
    fastq("in/barcode01/chunk_0.fastq.gz", [("a", "ACGT")])
    fastq("in/barcode01/chunk_1.fastq.gz", [("b", "ACGT")])
    fastq("in/unclassified/x.fastq.gz", [("c", "ACGT")])
    fasta("in/asm.fasta", [("c1", "ACGT")])
    (d / "README.md").write_text("ignored")
    samples = {s.name: s for s in discover(d)}
    assert sorted(samples) == ["asm", "barcode01"]
    assert len(samples["barcode01"].files) == 2
    assert samples["asm"].fmt == "fasta"


def test_duplicate_names_rejected(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/s1.fastq.gz", [("a", "ACGT")])
    fastq("in/s1.fq", [("a", "ACGT")], gz=False)
    with pytest.raises(RasherError, match="same sample name"):
        discover(tmp_path / "in")


def test_mixed_formats_rejected(tmp_path, fastq, fasta):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/a.fastq.gz", [("a", "ACGT")])
    fasta("in/s/b.fasta", [("b", "ACGT")])
    with pytest.raises(RasherError, match="mixes fasta and fastq"):
        discover(tmp_path / "in")


def test_empty_files_skipped_and_nothing_left_is_an_error(tmp_path):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "empty.fastq").write_text("")
    with pytest.raises(RasherError, match="No fasta/fastq"):
        discover(tmp_path / "in")


@pytest.mark.parametrize("name", ["bad name", "Reference", "all_assemblies", "s1\n"])
def test_invalid_or_reserved_names(tmp_path, fastq, name):
    fastq(f"{name}.fastq.gz", [("a", "ACGT")])
    with pytest.raises(RasherError, match="Invalid sample name|reserved"):
        discover(tmp_path / f"{name}.fastq.gz")


def test_sample_sheet(tmp_path, fastq):
    fastq("a1.fastq.gz", [("a", "ACGT")])
    fastq("a2.fastq.gz", [("b", "ACGT")])
    fastq("b.fastq.gz", [("c", "ACGT")])
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("# comment\nSample,File\nA,a1.fastq.gz;a2.fastq.gz\nB,b.fastq.gz\nA,a2.fastq.gz\n")
    samples = {s.name: s for s in read_sample_sheet(sheet)}
    assert [f.name for f in samples["A"].files] == ["a1.fastq.gz", "a2.fastq.gz"]  # a2 listed twice: used once
    assert samples["B"].files == [tmp_path / "b.fastq.gz"]


def test_sample_sheet_paths_keep_their_inner_spaces(tmp_path, fastq):
    (tmp_path / "reads").mkdir()
    fastq("reads/a  b.fastq.gz", [("a", "ACGT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\nab \t reads/a  b.fastq.gz \n")
    [sample] = read_sample_sheet(sheet)
    assert sample.name == "ab" and sample.files == [tmp_path / "reads" / "a  b.fastq.gz"]


def test_sample_sheet_errors(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("name\tpath\nA\tx.fq\n")
    with pytest.raises(RasherError, match="needs the columns"):
        read_sample_sheet(sheet)
    sheet.write_text("sample\tfile\nA\tmissing.fq\n")
    with pytest.raises(RasherError, match="file not found"):
        read_sample_sheet(sheet)


def test_files_in_hidden_subfolders_are_ignored(tmp_path, fastq):
    from rasher.samples import discover
    (tmp_path / "in" / "barcode01" / ".hidden").mkdir(parents=True)
    for path in ("in/barcode01/a.fastq", "in/barcode01/.hidden/b.fastq"):
        (tmp_path / path).write_text("@r\nACGT\n+\nIIII\n")
    [sample] = discover(tmp_path / "in")
    assert [f.name for f in sample.files] == ["a.fastq"]


def test_tsv_sheets_quotes_and_comments(tmp_path):
    from rasher.metadata import sheet_metadata
    for name in ("a.fastq", "a b.fastq"):
        (tmp_path / name).write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.tsv"
    # Cells entirely in quotes (a spreadsheet's TSV export) lose them, in the header too, and "" inside is a quote;
    # a line starting with # after the header is a comment (a sample left out), for the sheet's metadata too
    sheet.write_text('"sample"\t"file"\t"note"\n"a"\t"a.fastq"\t""\n"b"\t"a b.fastq"\t"5"" tube, ""fine"""\n'
                     "#c\ta.fastq\tx\n")
    samples = read_sample_sheet(sheet)
    assert [(s.name, [f.name for f in s.files]) for s in samples] == [("a", ["a.fastq"]), ("b", ["a b.fastq"])]
    assert sheet_metadata(sheet).rows == {"a": {"note": ""}, "b": {"note": '5" tube, "fine"'}}
    sheet.write_text("sample\tfile\na\ta.fastq\n# b\ta b.fastq\n")
    assert [s.name for s in read_sample_sheet(sheet)] == ["a"]
    sheet.write_text("sample\tfile\tnote\tnote\na\ta.fastq\tx\ty\n")  # Still an error, saying what to do
    with pytest.raises(RasherError, match="column name 'note' is used more than once .* give each column"):
        read_sample_sheet(sheet)


def test_csv_sheets_quotes_and_unclosed_quotes(tmp_path):
    for name in ("s1", "s2", "a b"):
        (tmp_path / f"{name}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.csv"
    # A space after a closing quote is kept, then stripped
    sheet.write_text('sample,file\n"s1" ,s1.fastq\ns2,"a b.fastq" \n')
    assert [(s.name, [f.name for f in s.files]) for s in read_sample_sheet(sheet)] == [
        ("s1", ["s1.fastq"]), ("s2", ["a b.fastq"])]
    sheet.write_text('sample,file\ns1,"s1.fastq\ns2,s2.fastq\n')  # Not s2 swallowed into s1's file name
    with pytest.raises(RasherError, match="line 2: an unclosed quote"):
        read_sample_sheet(sheet)


def test_sample_sheet_in_windows_1252_is_read_with_a_warning(tmp_path, caplog):
    (tmp_path / "\u00e9t\u00e9.fastq").write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.csv"
    sheet.write_bytes("sample,file\ns1,\u00e9t\u00e9.fastq\n".encode("cp1252"))
    assert [f.name for f in read_sample_sheet(sheet)[0].files] == ["\u00e9t\u00e9.fastq"]
    assert f"Sample sheet {sheet} is not UTF-8" in caplog.text


@pytest.mark.parametrize("name", [".", "..", ".x"])
def test_sample_names_cannot_start_with_a_dot(tmp_path, fastq, name):
    """'..' would name the output folder, '.' the assembly step folder: both are deleted before an assembly."""
    fastq("a.fastq.gz", [("a", "ACGT")])
    sheet = tmp_path / "s.tsv"
    sheet.write_text(f"sample\tfile\n{name}\ta.fastq.gz\n")
    with pytest.raises(RasherError, match=f"Invalid sample name '{name}' .*not starting with a dot"):
        read_sample_sheet(sheet)
    if name == ".x":  # A file in a folder with such a name is hidden (skipped); a single file is checked
        with pytest.raises(RasherError, match="Invalid sample name '.x'"):
            discover(fastq(".x.fastq.gz", [("a", "ACGT")]))


def test_linked_subfolders_are_followed_and_folders_without_reads_are_reported(tmp_path, fastq, caplog):
    import logging
    (tmp_path / "elsewhere" / "pass").mkdir(parents=True)
    fastq("elsewhere/pass/chunk_0.fastq.gz", [("a", "ACGT")])
    fastq("elsewhere/pass/chunk_1.fastq.gz", [("b", "ACGT")])
    d = tmp_path / "in"
    (d / "barcode01").mkdir(parents=True)
    (d / "barcode01" / "pass").symlink_to(tmp_path / "elsewhere" / "pass")  # A linked subfolder
    (d / "barcode01" / "again").symlink_to(tmp_path / "elsewhere" / "pass")  # The same folder twice: read once
    (d / "barcode01" / "loop").symlink_to(d / "barcode01")  # A loop
    fastq("in/barcode01/chunk_2.fastq.gz", [("c", "ACGT")])
    (d / "barcode02").mkdir()
    (d / "barcode02" / "notes.txt").write_text("no reads")
    (d / "barcode03").symlink_to(tmp_path / "elsewhere")  # A linked sample folder
    with caplog.at_level(logging.INFO, logger="rasher"):
        samples = {s.name: s for s in discover(d)}
    assert sorted(samples) == ["barcode01", "barcode03"]
    assert [f.name for f in samples["barcode01"].files] == ["chunk_0.fastq.gz", "chunk_1.fastq.gz",
                                                             "chunk_2.fastq.gz"]
    assert [f.name for f in samples["barcode03"].files] == ["chunk_0.fastq.gz", "chunk_1.fastq.gz"]
    assert f"Skipping folder {d / 'barcode02'}: no fasta/fastq files" in caplog.text


# ----- Read pairing ---------------------------------------------------------------------------------------------

def test_pair_stem_patterns():
    from rasher.samples import split_pair_stem
    assert split_pair_stem("x_R1.fastq.gz") == ("x", ("_", "R", ""), 1)
    assert split_pair_stem("x_r2.fq") == ("x", ("_", "R", ""), 2)
    assert split_pair_stem("x_1.fastq.gz") == ("x", ("_", "", ""), 1)
    assert split_pair_stem("x.R2.fasta") == ("x", (".", "R", ""), 2)
    assert split_pair_stem("x_S1_L001_R1_001.fastq.gz") == ("x_S1_L001", ("_", "R", "_001"), 1)
    assert split_pair_stem("a_1_2.fastq") == ("a_1", ("_", "", ""), 2)  # The last marker wins
    assert split_pair_stem("plain.fastq.gz") is None
    assert split_pair_stem("x_R3.fastq.gz") is None
    assert split_pair_stem("notaread.txt") is None


def test_flat_folder_pairs_one_sample(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/x_R2.fastq.gz", [("r", "TTTT")])
    fastq("in/y.fastq.gz", [("s", "ACGT")])
    samples = {s.name: s for s in discover(tmp_path / "in")}
    assert sorted(samples) == ["x", "y"]
    assert [(p1.name, p2.name) for p1, p2 in samples["x"].pairs] == [("x_R1.fastq.gz", "x_R2.fastq.gz")]
    assert [f.name for f in samples["x"].files] == ["x_R1.fastq.gz", "x_R2.fastq.gz"]
    assert samples["y"].pairs == []


def test_sra_style_underscore_numbers_pair(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/Acc_ERR123_1.fastq.gz", [("r", "ACGT")])
    fastq("in/Acc_ERR123_2.fastq.gz", [("r", "TTTT")])
    samples = discover(tmp_path / "in")
    assert [s.name for s in samples] == ["Acc_ERR123"] and samples[0].pairs


def test_r1_and_underscore_2_do_not_pair(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/x_2.fastq.gz", [("r", "TTTT")])
    samples = {s.name for s in discover(tmp_path / "in")}
    assert samples == {"x_R1", "x_2"}  # Different marker styles: two single-end samples


def test_file_without_a_mate_keeps_its_full_stem(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/sample_v2.fastq.gz", [("r", "ACGT")])
    assert [s.name for s in discover(tmp_path / "in")] == ["sample_v2"]


def test_duplicate_mates_stay_single(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/x_R1.fq", [("r", "ACGT")], gz=False)
    fastq("in/x_R2.fastq.gz", [("r", "TTTT")])
    with pytest.raises(RasherError, match="same sample name"):  # Both x_R1 files claim the same name
        discover(tmp_path / "in")


def test_subfolder_files_pair(tmp_path, fastq):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/runA_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/s/runA_R2.fastq.gz", [("r", "TTTT")])
    fastq("in/s/runB_R1.fastq.gz", [("q", "ACGT")])
    fastq("in/s/runB_R2.fastq.gz", [("q", "TTTT")])
    (sample,) = discover(tmp_path / "in")
    assert sample.name == "s" and len(sample.pairs) == 2 and len(sample.files) == 4


def test_orphan_next_to_pairs_is_skipped_with_a_warning(tmp_path, fastq, caplog):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/a_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/s/a_R2.fastq.gz", [("r", "TTTT")])
    fastq("in/s/orphan.fastq.gz", [("q", "AAAA")])
    (sample,) = discover(tmp_path / "in")
    assert len(sample.pairs) == 1 and [f.name for f in sample.files] == ["a_R1.fastq.gz", "a_R2.fastq.gz"]
    assert "orphan.fastq.gz" in caplog.text and "skipped" in caplog.text


def test_all_unpaired_subfolder_is_single_end(tmp_path, fastq):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/chunk_a.fastq.gz", [("r", "ACGT")])
    fastq("in/s/chunk_b.fastq.gz", [("q", "TTTT")])
    (sample,) = discover(tmp_path / "in")
    assert sample.pairs == [] and len(sample.files) == 2


def test_same_names_in_different_folders_do_not_cross_pair(tmp_path, fastq):
    (tmp_path / "in" / "s" / "runA").mkdir(parents=True)
    (tmp_path / "in" / "s" / "runB").mkdir(parents=True)
    fastq("in/s/runA/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/s/runB/x_R2.fastq.gz", [("r", "TTTT")])
    (sample,) = discover(tmp_path / "in")
    assert sample.pairs == []  # Mates pair within one folder only


def test_sheet_files_pair(tmp_path, fastq):
    fastq("x_1.fastq.gz", [("r", "ACGT")])
    fastq("x_2.fastq.gz", [("r", "TTTT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\np\tx_1.fastq.gz;x_2.fastq.gz\n")
    (sample,) = read_sample_sheet(sheet)
    assert sample.name == "p" and len(sample.pairs) == 1


def test_mates_that_do_not_pair_are_single_end_with_a_warning(tmp_path, fastq, caplog):
    (tmp_path / "in" / "dup").mkdir(parents=True)
    fastq("in/dup/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/dup/x_R1.fastq", [("r", "ACGT")], gz=False)  # A duplicate mate
    fastq("in/dup/x_R2.fastq.gz", [("r", "TTTT")])
    (sample,) = discover(tmp_path / "in")
    assert sample.pairs == [] and len(sample.files) == 3
    assert ("Sample dup: x_R1.fastq, x_R1.fastq.gz, x_R2.fastq.gz look like read mates but were not paired (a "
            "duplicate mate file, or mates in different folders): processed as single-end reads") in caplog.text
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    fastq("a/y_R1.fastq.gz", [("r", "ACGT")])
    fastq("b/y_R2.fastq.gz", [("r", "TTTT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\ny\ta/y_R1.fastq.gz;b/y_R2.fastq.gz\n")
    caplog.clear()
    (sample,) = read_sample_sheet(sheet)
    assert sample.pairs == [] and "Sample y: y_R1.fastq.gz, y_R2.fastq.gz look like read mates" in caplog.text


def test_fasterq_dump_unpaired_reads_next_to_their_pair_are_skipped(tmp_path, fastq, caplog):
    (tmp_path / "in").mkdir()
    fastq("in/SRR1.fastq", [("u", "ACGT")], gz=False)  # fasterq-dump --split-3
    fastq("in/SRR1_1.fastq", [("r", "ACGT")], gz=False)
    fastq("in/SRR1_2.fastq", [("r", "TTTT")], gz=False)
    fastq("in/SRR2.fastq", [("s", "ACGT")], gz=False)
    samples = {s.name: s for s in discover(tmp_path / "in")}
    assert sorted(samples) == ["SRR1", "SRR2"] and len(samples["SRR1"].pairs) == 1 and samples["SRR2"].files
    assert (f"Skipping {tmp_path / 'in' / 'SRR1.fastq'}: the sample SRR1 is the pair SRR1_1.fastq + SRR1_2.fastq"
            in caplog.text)


def test_a_file_does_not_pair_with_a_link_to_itself(tmp_path, fastq, caplog):
    d = tmp_path / "in"
    d.mkdir()
    fastq("in/x_R1.fastq.gz", [("r", "ACGT")])
    (d / "x_R2.fastq.gz").symlink_to(d / "x_R1.fastq.gz")
    assert [(s.name, s.pairs) for s in discover(d)] == [("x_R1", []), ("x_R2", [])]
    assert (f"{d / 'x_R1.fastq.gz'} and {d / 'x_R2.fastq.gz'} are the same file (one links to the other): not "
            "paired") in caplog.text
    (d / "s").mkdir()
    fastq("in/s/y_R1.fastq.gz", [("r", "ACGT")])
    (tmp_path / "in" / "s" / "y_R2.fastq.gz").symlink_to("y_R1.fastq.gz")
    caplog.clear()
    sample = {s.name: s for s in discover(tmp_path / "in")}["s"]
    assert sample.pairs == [] and [f.name for f in sample.files] == ["y_R1.fastq.gz"]  # Its reads used once
    assert "y_R2.fastq.gz is the file" in caplog.text and "look like read mates" not in caplog.text


def test_mates_of_different_formats_do_not_pair(tmp_path, fastq, fasta):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fasta("in/x_R1.fasta", [("r", "ACGT")])
    fastq("in/x_R2.fastq", [("r", "TTTT")], gz=False)
    assert [(s.name, s.fmt, s.pairs) for s in discover(tmp_path / "in")] == [("x_R1", "fasta", []),
                                                                            ("x_R2", "fastq", [])]
    fasta("in/s/y_R1.fasta", [("r", "ACGT")])
    fastq("in/s/y_R2.fastq", [("r", "TTTT")], gz=False)
    with pytest.raises(RasherError, match="Sample s mixes fasta and fastq"):  # One sample: still an error
        discover(tmp_path / "in")


def test_a_single_mate_file_gets_a_pairing_hint(fastq, caplog):
    import logging
    with caplog.at_level(logging.INFO, logger="rasher"):
        (sample,) = discover(fastq("x_R1.fastq.gz", [("r", "ACGT")]))
    assert sample.name == "x_R1" and sample.pairs == []
    assert "x_R1.fastq.gz looks like one mate of a pair: give the folder or a sample sheet" in caplog.text


def test_single_end_sample_collides_with_samples_named_as_its_mates(tmp_path, fastq):
    """A single-end sample deletes the X_R1 + X_R2 reads of an earlier paired run of it."""
    for name in ("x", "x_R1", "x_R2", "y", "y_R1"):
        fastq(f"{name}.fastq.gz", [("r", "ACGT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\ny\ty.fastq.gz\ny_R1\ty_R1.fastq.gz\n")
    assert len(read_sample_sheet(sheet)) == 2  # One of the two: no X_R2 to delete with it
    sheet.write_text("sample\tfile\nx\tx.fastq.gz\nx_R1\tx_R1.fastq.gz\nx_R2\tx_R2.fastq.gz\n")
    with pytest.raises(RasherError, match="Sample names 'x', 'x_R1' and 'x_R2' collide"):
        read_sample_sheet(sheet)


def test_paired_sample_name_collides_with_mate_named_sample(tmp_path, fastq):
    (tmp_path / "in" / "x_R1").mkdir(parents=True)
    fastq("in/x_R1/reads.fastq.gz", [("q", "AAAA")])
    fastq("in/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/x_R2.fastq.gz", [("r", "TTTT")])
    with pytest.raises(RasherError, match="collide"):
        discover(tmp_path / "in")


def test_sample_names_that_differ_only_by_case_are_rejected(tmp_path, fastq):
    """On a case-insensitive file system (macOS), S1 and s1 would write the same files."""
    for name in ("a", "b", "x_R1", "x_R2"):
        fastq(f"{name}.fastq.gz", [("r", "ACGT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\nS1\ta.fastq.gz\ns1\tb.fastq.gz\n")
    with pytest.raises(RasherError, match="Sample names 'S1' .* and 's1' .* differ only by case"):
        read_sample_sheet(sheet)
    (tmp_path / "in").mkdir()
    fastq("in/S1.fastq.gz", [("r", "ACGT")])
    fastq("in/s1.fastq.gz", [("r", "ACGT")])
    if len(list((tmp_path / "in").iterdir())) == 2:  # On a case-insensitive file system they are one file
        with pytest.raises(RasherError, match="differ only by case"):
            discover(tmp_path / "in")
    # A paired sample x writes x_R1 files: those of a sample X_r1 on such a file system
    sheet.write_text("sample\tfile\nx\tx_R1.fastq.gz;x_R2.fastq.gz\nX_r1\ta.fastq.gz\n")
    with pytest.raises(RasherError, match="Sample names 'x' \\(paired\\) and 'X_r1' collide"):
        read_sample_sheet(sheet)


def test_a_leftover_of_the_other_format_next_to_pairs_is_skipped(tmp_path, fastq, fasta, caplog):
    """A stray fasta (contigs of an earlier assembly) next to fastq pairs is a leftover, not a mixed sample."""
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/x_R1.fastq.gz", [("r", "ACGT")])
    fastq("in/s/x_R2.fastq.gz", [("r", "TTTT")])
    fasta("in/s/contigs.fasta", [("c1", "ACGT")])
    (sample,) = discover(tmp_path / "in")
    assert sample.fmt == "fastq" and [f.name for f in sample.files] == ["x_R1.fastq.gz", "x_R2.fastq.gz"]
    assert "Sample s: contigs.fasta has no mate next to the 1 pair(s): skipped" in caplog.text
    # Pairs of both formats are one sample's reads: still an error
    fasta("in/s/y_R1.fasta", [("r", "ACGT")])
    fasta("in/s/y_R2.fasta", [("r", "TTTT")])
    with pytest.raises(RasherError, match="Sample s mixes fasta and fastq"):
        discover(tmp_path / "in")
