import pytest

from scram2Plot.utils import ScramCSV

CSV = (
    "Header,len,sRNA,Position,Strand,Times aligned,rep1,rep2\n"
    "chr1,100,AAAAA,10,+,1,4.0,6.0\n"
    "chr1,100,CCCCC,20,-,2,1.0,1.0\n"
    "chr2,200,GGGGG,30,+,1,10.0,10.0\n"
)


@pytest.fixture
def scram_csv(tmp_path):
    file_path = tmp_path / "aligned_21.csv"
    file_path.write_text(CSV)
    sc = ScramCSV()
    sc.load_csv(str(file_path))
    return sc


def test_load_csv_srna_len(scram_csv):
    assert scram_csv.srna_len == 21
    assert len(scram_csv.df) == 3


@pytest.mark.parametrize("header", ["001", "NA", "N/A", "null", "1e3"])
def test_load_csv_preserves_literal_headers(tmp_path, header):
    path = tmp_path / "literal_21.csv"
    path.write_text(
        "Header,len,sRNA,Position,Strand,Times aligned,rep1\n"
        f"{header},100,{'A' * 21},1,+,1,10\n"
    )
    sc = ScramCSV()
    sc.load_csv(str(path))
    assert sc.subset_data(header=header).df["Header"].tolist() == [header]


def test_subset_by_header_only(scram_csv):
    # used to raise KeyError dropping the never-created Mean_Count column
    sub = scram_csv.subset_data(header="chr1")
    assert len(sub.df) == 2
    assert "Mean_Count" not in sub.df.columns


def test_subset_chained(scram_csv):
    # filtering a filtered instance used to break on index misalignment
    sub = scram_csv.subset_data(header="chr1")
    sub2 = sub.subset_data(strand="+")
    assert len(sub2.df) == 1
    assert sub2.df.iloc[0]["sRNA"] == "AAAAA"


def test_subset_min_count(scram_csv):
    sub = scram_csv.subset_data(min_count=2.0)
    # mean counts: 5.0, 1.0, 10.0 -> keeps rows 1 and 3
    assert len(sub.df) == 2
    assert "Mean_Count" not in sub.df.columns


def test_subset_preserves_srna_len(scram_csv):
    sub = scram_csv.subset_data(header="chr2")
    assert sub.srna_len == 21
