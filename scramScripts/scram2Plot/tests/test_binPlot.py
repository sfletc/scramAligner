import numpy as np

from scram2Plot.binPlot import calculate_means_and_ses, create_dataframe, get_bin_data
from scram2Plot.profilePlot import SingleAlignment, SingleRefProfile


def make_profile(alignments, ref_len, replicates):
    profile = SingleRefProfile()
    profile.ref_len = ref_len
    profile.replicates = replicates
    profile.all_alignments = alignments
    return profile


def test_get_bin_data_bin_boundaries():
    bin_size = 100
    profile = make_profile(
        [
            SingleAlignment("A" * 21, 1, "+", 1, np.array([2.0])),  # first base -> bin 0
            SingleAlignment("A" * 21, 100, "+", 1, np.array([3.0])),  # last base of bin 0
            SingleAlignment("A" * 21, 101, "+", 1, np.array([5.0])),  # first base of bin 1
            SingleAlignment("A" * 21, 150, "-", 1, np.array([7.0])),  # reverse strand, bin 1
        ],
        ref_len=250,
        replicates=1,
    )
    bins = get_bin_data(profile, bin_size)
    assert bins.shape == (3, 2, 1)
    assert bins[0, 0, 0] == 5.0  # positions 1 and 100 both in bin 0
    assert bins[1, 0, 0] == 5.0
    assert bins[1, 1, 0] == 7.0
    assert bins[2].sum() == 0.0


def test_get_bin_data_empty_profile():
    profile = SingleRefProfile()
    profile.ref_len = 100
    profile.replicates = 2
    bins = get_bin_data(profile, 50)
    assert bins.shape == (2, 2, 2)
    assert bins.sum() == 0.0


def test_calculate_means_and_ses_sample_sd():
    bins = np.zeros((1, 2, 3))
    bins[0, 0] = [1.0, 2.0, 3.0]
    means, ses = calculate_means_and_ses(bins)
    assert means[0, 0] == 2.0
    expected_se = np.std([1.0, 2.0, 3.0], ddof=1) / np.sqrt(3)
    assert np.isclose(ses[0, 0], expected_se)


def test_calculate_means_and_ses_single_replicate():
    bins = np.ones((2, 2, 1))
    means, ses = calculate_means_and_ses(bins)
    assert not np.isnan(ses).any()
    assert (ses == 0).all()


def test_create_dataframe_negates_reverse_strand():
    bins = np.zeros((2, 2, 2))
    means = np.array([[1.0, 2.0], [3.0, 4.0]])
    ses = np.array([[0.1, 0.2], [0.3, 0.4]])
    df = create_dataframe(bins, means, ses)
    fwd = df[df["list"] == "first"]
    rvs = df[df["list"] == "second"]
    assert fwd["mean"].tolist() == [1.0, 3.0]
    assert rvs["mean"].tolist() == [-2.0, -4.0]
    assert rvs["se"].tolist() == [0.2, 0.4]
