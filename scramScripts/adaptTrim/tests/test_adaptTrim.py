import pytest
import src.adaptTrim as flex_adapter

def test_trim_adapter_fastq():
    """Test the trim_adapter function."""
    # Test with a read that has a perfect match to the adapter
    read = ["@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","TGTGACGACTCTCGGCAACGGATATCCTCATGGAATTCTCGGGTGCCAAGG","+","F,FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"]
    adapter = "TGGAATTCTCGGGTGCCAAGG"
    min_overlap = 8
    min_len = 18
    trim5 = 4
    trim3 = 4
    error=0.1
    trimmed_read = flex_adapter.trim_adapter(read, adapter, min_len, trim5, trim3, min_overlap, error)
    assert trimmed_read == [("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC","+","FFFFFFFFFFF:FFFFFFFFFF")]

# test loading a fastq file
def test_load_fastq():
    """Test the load_fastq function."""
    # Test with a fastq file with 4 reads
    fq = "tests/test_data/test.fq"
    loaded_fq = flex_adapter.load_fastq(fq)
    assert len(loaded_fq) == 80

# test rejecting a read with low quality
def test_reject_low_quality():
    """Test the reject_low_quality function."""
    # Test with a read that has a mean quality score above the error threshold
    trimmed_reads = [("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC","+","***************")]
    max_error = 0.1
    trimmed_reads = flex_adapter.reject_low_quality(trimmed_reads, max_error)
    assert len(trimmed_reads) == 0

#test retaining a read with high quality
def test_retain_high_quality():
    """Test the reject_low_quality function."""
    # Test with a read that has a mean quality score below the error threshold
    trimmed_reads = [("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC","+","***************")]
    max_error = 0.2
    trimmed_reads = flex_adapter.reject_low_quality(trimmed_reads, max_error)
    assert len(trimmed_reads) == 1


# test writing a fasta file
def test_write_fasta():
    """Test the write_fasta function."""
    # Test with a list of 2 reads
    trimmed_reads = [("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC"),("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC")]
    out_fa = "tests/test_data/test_write.fa"
    flex_adapter.write_fastx(trimmed_reads, out_fa, fa=True, gz=False)
    with open(out_fa, 'r') as f:
        assert len(f.readlines()) == 4

def test_write_fastq():
    """Test the write_fastq function."""
    # Test with a list of 2 reads
    trimmed_reads = [("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC","+","FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"),("@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","ACGACTCTCGGCAACGGATATC","+","FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF")]
    out_fq = "tests/test_data/test_write.fq"
    flex_adapter.write_fastx(trimmed_reads, out_fq, fa=False, gz=False)
    with open(out_fq, 'r') as f:
        assert len(f.readlines()) == 8

# test the CLI
def test_flex_adapter_trim():
    """Test the flex_adapter_trim function."""
    # Test with a fastq file with 4 reads
    fq = "tests/test_data/test.fq"
    adapter = "TGGAATTCTCGGGTGCCAAGG"
    out_fa = "tests/test_data/test.fa"
    flex_adapter.flex_adapter_trim(fq, adapter, out_fa)
    with open(out_fa, 'r') as f:
        assert len(f.readlines()) == 20

# test read that is too short

def test_trim_adapter_short_read():
    """Test the trim_adapter function with a read that is too short."""
    # Test with a read that has a perfect match to the adapter
    read = ["@A00599:365:HFVTWDSX5:1:1101:4435:1000 1:N:0:AGTCAAAT+GGGGGGGG","TGTTGGAATTCTCGGGTGCCAAGGGGCCTCATGGAATTCTCGGGTGCCAAGG","+","F,FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"]
    adapter = "TGGAATTCTCGGGTGCCAAGG"
    min_overlap = 8
    min_len = 18
    trim5 = 4
    trim3 = 4
    fa = True
    trimmed_read = flex_adapter.trim_adapter(read, adapter, min_len, trim5, trim3, min_overlap, fa)
    assert trimmed_read == []