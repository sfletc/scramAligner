import numpy as np
import pytest

from scram2Plot.profilePlot import (
    DNA,
    AlignmentPlot,
    DataForPlot,
    RefProfiles,
    SingleAlignment,
    SingleRefProfile,
    color_for_length,
)


def test_init():
    dna = DNA("AGCT")
    assert dna.sequence == "AGCT"


def test_init_uppercase():
    dna = DNA("agct")
    assert dna.sequence == "AGCT"


def test_len():
    dna = DNA("AGCT")
    assert len(dna) == 4


def test_getitem():
    dna = DNA("AGCT")
    assert dna[0] == "A"
    assert dna[1] == "G"
    assert dna[2] == "C"
    assert dna[3] == "T"


def test_hash():
    dna1 = DNA("AGCT")
    dna2 = DNA("AGCT")
    assert hash(dna1) == hash(dna2)


def test_repr():
    dna = DNA("AGCT")
    assert repr(dna) == "AGCT"


def test_eq():
    dna1 = DNA("AGCT")
    dna2 = DNA("AGCT")
    dna3 = DNA("AGCA")
    assert dna1 == dna2
    assert dna1 != dna3


def test_single_alignment_init():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.srna == "ACGT"
    assert s.position == 10
    assert s.strand == "+"
    assert s.times_aligned == 5
    assert np.array_equal(s.indv_alignments, [1, 2, 3, 4, 5])


def test_single_alignment_srna_len():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.srna_len() == 4


def test_single_alignment_mean_alignments():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s.mean_alignments() == 3.0


def test_single_alignment_str():
    s = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    assert str(s) == "ACGT\t10\t+\t5\t[1, 2, 3, 4, 5]"


def test_single_alignment_eq():
    s1 = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    s2 = SingleAlignment("ACGT", 10, "+", 5, [1, 2, 3, 4, 5])
    s3 = SingleAlignment("ACGT", 10, "-", 5, [1, 2, 3, 4, 5])
    s4 = SingleAlignment("ACGT", 10, "+", 6, [1, 2, 3, 4, 5, 6])
    s5 = SingleAlignment("TGCA", 10, "+", 5, [1, 2, 3, 4, 5])
    assert s1 == s2
    assert s1 != s3
    assert s1 != s4
    assert s1 != s5


def test_single_ref_profile_init():
    s = SingleRefProfile()
    assert s.ref_len == 0
    assert s.all_alignments == []


def test_single_ref_profile_str():
    s = SingleRefProfile()
    assert str(s) == "0\t[]"


def test_single_ref_profile_eq():
    s1 = SingleRefProfile()
    s2 = SingleRefProfile()
    s3 = SingleRefProfile()
    s3.ref_len = 10
    s4 = SingleRefProfile()
    s4.all_alignments = [1, 2, 3]
    assert s1 == s2
    assert s1 != s3
    assert s1 != s4


def test_ref_profiles_init():
    r = RefProfiles()
    assert r.srna_len == 0
    assert r.replicates == 0
    assert r.single_ref_profiles == {}


def test_ref_profiles_eq():
    r1 = RefProfiles()
    r2 = RefProfiles()
    r3 = RefProfiles()
    r3.srna_len = 20
    r4 = RefProfiles()
    r4.replicates = 5
    r5 = RefProfiles()
    r5.single_ref_profiles["ref1"] = SingleRefProfile()
    r5.single_ref_profiles["ref2"] = SingleRefProfile()
    assert r1 == r2
    assert r1 != r3
    assert r1 != r4
    assert r1 != r5


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

CSV_HEADER = "Header,len,sRNA,Position,Strand,Times aligned,rep1,rep2,rep3\n"


def write_csv(tmp_path, rows, header=CSV_HEADER, name="profile_4.csv"):
    file_path = tmp_path / name
    file_path.write_text(header + "".join(rows))
    return str(file_path)


def test_ref_profiles_load_single_ref_profiles(tmp_path):
    r = RefProfiles()
    file_path = write_csv(
        tmp_path,
        [
            "r1,20,ACGT,10,+,5,1.0,2.0,3.0\n",
            "r1,20,ACGT,20,-,5,4.0,5.0,6.0\n",
            "r2,25,CGTA,5,+,5,7.0,8.0,9.0\n",
        ],
    )
    r.load_single_ref_profiles(file_path)
    assert r.srna_len == 4
    assert r.replicates == 3
    assert "r1" in r.single_ref_profiles
    assert "r2" in r.single_ref_profiles
    assert r.single_ref_profiles["r1"].ref_len == 20
    assert len(r.single_ref_profiles["r1"].all_alignments) == 2
    assert r.single_ref_profiles["r1"].all_alignments[0] == SingleAlignment(
        DNA("ACGT"), 10, "+", 5, [1.0, 2.0, 3.0]
    )
    assert r.single_ref_profiles["r1"].all_alignments[1] == SingleAlignment(
        DNA("ACGT"), 20, "-", 5, [4.0, 5.0, 6.0]
    )
    assert r.single_ref_profiles["r2"].ref_len == 25
    assert len(r.single_ref_profiles["r2"].all_alignments) == 1
    assert r.single_ref_profiles["r2"].all_alignments[0] == SingleAlignment(
        DNA("CGTA"), 5, "+", 5, [7.0, 8.0, 9.0]
    )


def test_load_rejects_mean_se_format(tmp_path):
    r = RefProfiles()
    file_path = write_csv(
        tmp_path,
        ["r1,20,ACGT,10,+,0.5,0.1,5\n"],
        header="Header,len,sRNA,Position,Strand,Count,Std. Err,Times aligned\n",
    )
    with pytest.raises(ValueError, match="mean/standard-error"):
        r.load_single_ref_profiles(file_path)


def test_load_missing_header_raises(tmp_path):
    r = RefProfiles()
    file_path = write_csv(tmp_path, ["r1,20,ACGT,10,+,5,1.0,2.0,3.0\n"])
    with pytest.raises(ValueError, match="no alignments"):
        r.load_single_ref_profiles(file_path, header="no_such_header")


@pytest.mark.parametrize("header", ["001", "NA", "N/A", "null", "1e3"])
def test_load_preserves_literal_headers(tmp_path, header):
    path = write_csv(tmp_path, [f"{header},100,AAAA,10,+,1,1.0,2.0,3.0\n"])
    profiles = RefProfiles()
    profiles.load_single_ref_profiles(path, header=header)
    assert list(profiles.single_ref_profiles) == [header]


@pytest.mark.parametrize("missing", ["", "NA", "NaN"])
def test_load_still_rejects_missing_counts(tmp_path, missing):
    path = write_csv(tmp_path, [f"NA,100,AAAA,10,+,1,1.0,{missing},3.0\n"])
    with pytest.raises(ValueError, match="missing or non-numeric"):
        RefProfiles().load_single_ref_profiles(path)


def test_load_subset_keeps_left_edge_overlap(tmp_path):
    # A read starting before the window but overlapping it must be kept, with a
    # negative relative position.
    r = RefProfiles()
    file_path = write_csv(
        tmp_path,
        [
            "r1,100,AAAA,48,+,1,1.0,1.0,1.0\n",  # overlaps window start at 50
            "r1,100,CCCC,10,+,1,1.0,1.0,1.0\n",  # fully before the window
            "r1,100,GGGG,60,+,1,1.0,1.0,1.0\n",  # inside
            "r1,100,TTTT,90,+,1,1.0,1.0,1.0\n",  # after the window
        ],
    )
    r.load_single_ref_profiles(file_path, header="r1", start=50, end=80)
    profile = r.single_ref_profiles["r1"]
    assert sorted(profile.srnas) == ["AAAA", "GGGG"]
    assert sorted(profile.positions.tolist()) == [-2, 10]


# ---------------------------------------------------------------------------
# DataForPlot transforms
# ---------------------------------------------------------------------------


def make_spd(arr_rows, srna_len, replicates, margin=0):
    spd = object.__new__(DataForPlot)
    spd.srna_len = srna_len
    spd.replicates = replicates
    spd.ref_len = arr_rows - 1
    spd.margin = margin
    spd.fwd = np.zeros((margin + arr_rows, replicates))
    spd.rvs = np.zeros((margin + arr_rows, replicates))
    spd.x_axis = np.arange(arr_rows)
    spd.y_flat = []
    return spd


def test_coverage_sum_overlapping_reads():
    # two overlapping reads: coverage must sum where they overlap
    spd = make_spd(20, 5, 1)
    arr = np.zeros((20, 1))
    arr[2] = 10.0
    arr[4] = 3.0
    cov = spd._coverage_per_strand(arr, abund=False)
    expected = np.zeros((20, 1))
    expected[2:7] += 10.0
    expected[4:9] += 3.0
    assert np.array_equal(cov, expected)
    # total coverage mass equals sum(counts) * srna_len
    assert cov.sum() == (10.0 + 3.0) * 5


def test_coverage_abund_winner_takes_all():
    # overlap resolved by the higher-mean read from either side
    spd = make_spd(20, 5, 1)
    arr = np.zeros((20, 1))
    arr[2] = 10.0  # high-abundance read at 2-6
    arr[4] = 3.0  # low-abundance read at 4-8
    cov = spd._coverage_per_strand(arr, abund=True)
    expected = np.zeros((20, 1))
    expected[2:7] = 10.0  # winner over the whole overlap
    expected[7:9] = 3.0  # remainder of the low read
    assert np.array_equal(cov, expected)

    # mirrored: low read first, high read second
    arr2 = np.zeros((20, 1))
    arr2[2] = 3.0
    arr2[4] = 10.0
    cov2 = spd._coverage_per_strand(arr2, abund=True)
    expected2 = np.zeros((20, 1))
    expected2[2:4] = 3.0
    expected2[4:9] = 10.0
    assert np.array_equal(cov2, expected2)


def test_coverage_short_reference():
    # reference shorter than the read length: no valid windows, all zeros
    spd = make_spd(3, 5, 2)
    arr = np.ones((3, 2))
    cov = spd._coverage_per_strand(arr, abund=False)
    assert np.array_equal(cov, np.zeros((3, 2)))


def test_error_bounds_uses_sample_sd():
    spd = make_spd(4, 21, 3)
    arr = np.array([[1.0, 2.0, 3.0]] * 4)
    bounds = spd._error_bounds(arr)
    se = np.std([1.0, 2.0, 3.0], ddof=1) / np.sqrt(3)
    assert np.allclose(bounds[:, 0], 2.0 - se)
    assert np.allclose(bounds[:, 1], 2.0 + se)


def test_error_bounds_single_replicate():
    spd = make_spd(4, 21, 1)
    arr = np.array([[5.0]] * 4)
    bounds = spd._error_bounds(arr)
    assert not np.isnan(bounds).any()
    assert np.array_equal(bounds[:, 0], bounds[:, 1])


def test_flatten_block_order_and_idempotent():
    spd = make_spd(10, 21, 3)
    rng = np.random.default_rng(1)
    spd.fwd = rng.random((10, 3))
    spd.rvs = rng.random((10, 3))
    spd.flatten()
    spd.flatten()  # must not accumulate
    assert len(spd.y_flat) == 6
    for i in range(3):
        # forward columns first, matched (negated) reverse columns after
        assert np.array_equal(spd.y_flat[i], spd.fwd[:, i])
        assert np.array_equal(spd.y_flat[i + 3], -spd.rvs[:, i])


def test_smoothing_keeps_peak_centred():
    data = np.zeros(101)
    data[50] = 100.0
    smoothed = DataForPlot._smoothTriangle(data, 5)
    assert np.argmax(smoothed) == 50
    assert np.isclose(smoothed.sum(), 100.0)


@pytest.mark.parametrize("length, degree", [(1, 5), (5, 4), (8, 5), (11, 5)])
def test_smoothing_short_data_preserves_length_and_values(length, degree):
    data = np.full(length, 7.0)
    smoothed = DataForPlot._smoothTriangle(data, degree)
    assert smoothed.shape == data.shape
    assert np.allclose(smoothed, data)


def test_smoothing_short_plot_saves(tmp_path, monkeypatch):
    monkeypatch.setattr("matplotlib.pyplot.show", lambda: None)
    write_csv(tmp_path, ["r1,8,AAAA,1,+,1,1.0,2.0,3.0\n"], name="short_4.csv")
    AlignmentPlot(
        str(tmp_path / "short"), [4], "r1", smoothing_window=5,
        cov=True, abund=False, se=False,
    ).generate_plot()
    assert (tmp_path / "short_r1.svg").stat().st_size > 0


def test_data_for_plot_round_trip(tmp_path):
    # synthetic CSV -> DataForPlot -> coverage; total mass must be conserved
    r = RefProfiles()
    file_path = write_csv(
        tmp_path,
        [
            "r1,50,AAAAA,3,+,1,2.0,4.0,6.0\n",
            "r1,50,CCCCC,5,+,1,1.0,1.0,1.0\n",
            "r1,50,GGGGG,10,-,1,3.0,3.0,3.0\n",
        ],
        name="profile_5.csv",
    )
    r.load_single_ref_profiles(file_path)
    spd = DataForPlot(r, "r1")
    fwd_mass = spd.fwd.sum()
    rvs_mass = spd.rvs.sum()
    assert fwd_mass == 2.0 + 4.0 + 6.0 + 3.0
    assert rvs_mass == 9.0
    spd.convert_to_coverage(abund=False)
    assert spd.fwd.sum() == fwd_mass * 5
    assert spd.rvs.sum() == rvs_mass * 5
    spd.convert_to_error_bounds()
    spd.flatten()
    assert len(spd.y_flat) == 4  # lower/upper for each strand
    assert all(len(col) == spd.ref_len + 1 for col in spd.y_flat)


def test_subset_left_edge_coverage(tmp_path):
    # Coverage inside a subset window must match the same positions from a
    # whole-reference run, including bases covered by a read that starts
    # before the window.
    rows = [
        "r1,100,AAAAA,48,+,1,1.0\n",  # covers 48-52, overlaps window start 50
        "r1,100,CCCCC,60,+,1,2.0\n",
    ]
    header_1rep = "Header,len,sRNA,Position,Strand,Times aligned,rep1\n"

    r_full = RefProfiles()
    r_full.load_single_ref_profiles(
        write_csv(rows=rows, tmp_path=tmp_path, header=header_1rep, name="full_5.csv")
    )
    spd_full = DataForPlot(r_full, "r1")
    spd_full.convert_to_coverage(abund=False)

    start, end = 50, 80
    r_sub = RefProfiles()
    r_sub.load_single_ref_profiles(
        write_csv(rows=rows, tmp_path=tmp_path, header=header_1rep, name="sub_5.csv"),
        header="r1",
        start=start,
        end=end,
    )
    spd_sub = DataForPlot(r_sub, "r1")
    spd_sub.convert_to_coverage(abund=False)

    # subset index i corresponds to full-reference index start + i
    assert np.array_equal(
        spd_sub.fwd[: end - start + 1], spd_full.fwd[start : end + 1]
    )


@pytest.mark.parametrize("abund", [False, True])
@pytest.mark.parametrize("start,end,padding", [(50, 100, 0), (50, 55, 0), (90, 90, 2)])
def test_subset_coverage_matches_full_at_both_edges(tmp_path, abund, start, end, padding):
    rows = [
        f"r1,200,{'A' * 21},{pos},{strand},1,{i + 1},{i + 2},{i + 3}\n"
        for i, pos in enumerate([29, 31, 48, 50, 60, 80, 90, 100, 101])
        for strand in ["+", "-"]
    ]
    path = write_csv(tmp_path, rows, name="edges_21.csv")
    full = RefProfiles()
    full.load_single_ref_profiles(path)
    full_data = DataForPlot(full, "r1", abund=abund)
    full_data.convert_to_coverage(abund=abund)

    subset = RefProfiles()
    subset.load_single_ref_profiles(path, header="r1", start=start, end=end, padding=padding)
    subset_data = DataForPlot(subset, "r1", abund=abund)
    subset_data.convert_to_coverage(abund=abund)
    subset_data.flatten()
    for strand in ["fwd", "rvs"]:
        expected = getattr(full_data, strand)[start - padding : end + padding + 1]
        assert np.array_equal(getattr(subset_data, strand), expected)
    assert all(len(y) == len(subset_data.x_axis) for y in subset_data.y_flat)


def test_subset_without_coverage_trims_both_margins(tmp_path):
    path = write_csv(tmp_path, [f"r1,200,{'A' * 21},100,+,1,1,2,3\n"])
    profiles = RefProfiles()
    profiles.load_single_ref_profiles(path, start=50, end=100)
    data = DataForPlot(profiles, "r1")
    data.trim_margin()
    data.flatten()
    assert data.fwd.shape == (51, 3)
    assert data.fwd[-1].tolist() == [1, 2, 3]
    assert all(len(y) == len(data.x_axis) for y in data.y_flat)


def test_color_for_length_fallback():
    assert color_for_length(21) == "red"
    # lengths outside the map must not raise
    assert color_for_length(35) is not None
    assert color_for_length(15) is not None
