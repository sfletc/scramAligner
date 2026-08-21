import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import argparse
import sys

# Shared color scheme for sRNA lengths; lengths outside the map get a fallback
# colormap color instead of raising KeyError.
SRNA_LEN_COLS = {
    18: "#f781bf",
    19: "#a65628",
    20: "#984ea3",
    21: "red",
    22: "blue",
    23: "#999999",
    24: "darkgreen",
    25: "brown",
    26: "#dede00",
    27: "orange",
    28: "yellow",
}


def color_for_length(srna_len):
    """Returns the plot color for an sRNA length, falling back to a colormap
    color for lengths outside the standard 18-28 nt range."""
    try:
        return SRNA_LEN_COLS[int(srna_len)]
    except (KeyError, ValueError):
        cmap = plt.get_cmap("tab10")
        return cmap(hash(str(srna_len)) % 10)


class DNA(object):
    """
    The DNA class represents a DNA sequence in a biological context.

    Class Attributes:
    dna_alphabet : set
        The nucleotides that can be found in a DNA sequence (Adenine, Guanine, Cytosine, Thymine, and N for unknown).

    Instance Attributes:
    sequence : str
        A string representing the DNA sequence. All characters are automatically converted to upper-case.

    Methods:
    __len__:
        Returns the length of the DNA sequence.
    __getitem__:
        Allows indexing into the DNA sequence.
    __hash__:
        Returns a hash value of the DNA sequence.
    __repr__:
        Returns a string representation of the DNA sequence.
    __eq__:
        Determines if two DNA objects have the same DNA sequence.

    Note: The DNA sequence should only contain characters that are in the dna_alphabet.
    No checks for invalid characters are performed.
    """

    dna_alphabet = set("AGCTN")

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence


class SingleAlignment(object):
    """
    The SingleAlignment class represents a single alignment of a small RNA (sRNA) sequence against a reference sequence.

    Each instance of the class reflects an exact match and carries a list of floating-point counts representing the number of sRNAs
    from each input read file that align at this location.

    Attributes:
    srna : DNA
        The DNA object representing the sRNA read.
    position : int
        The position in the reference sequence where the sRNA read aligns.
    strand : str
        The DNA strand ('+' or '-') where the sRNA read aligns.
    times_aligned : int
        The total number of positions the sRNA read can align to across all reference sequences.
    indv_alignments : list of float
        A numpy array of floating-point counts representing sRNA count for each individual sample.  Float as may be reads per million reads

    Methods:
    srna_len():
        Returns the length of the sRNA read.
    standard_error():
        Calculates the standard error of the individual alignment counts.
    mean_alignments():
        Calculates the mean of the individual alignment counts
    __str__():
        Returns a string representation of the SingleAlignment object, suitable for output.
    __repr__():
        Returns the same as __str__(), a string representation of the object.
    __eq__():
        Determines if two SingleAlignment objects are identical, including their individual alignments.
    """

    def __init__(self, srna, position, strand, times_aligned, indv_alignments):
        self.srna = srna
        self.position = position
        self.strand = strand
        self.times_aligned = times_aligned
        self.indv_alignments = indv_alignments

    def srna_len(self):
        return len(self.srna)

    def standard_error(self):  # TODO: check - likely don't neeed here
        return np.std(self.indv_alignments, ddof=1) / np.sqrt(
            np.size(self.indv_alignments)
        )

    def mean_alignments(self):
        return np.mean(self.indv_alignments)

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.srna,
            self.position,
            self.strand,
            self.times_aligned,
            self.indv_alignments,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.srna == other.srna
            and self.position == other.position
            and self.strand == other.strand
            and self.times_aligned == other.times_aligned
            and np.array_equal(self.indv_alignments, other.indv_alignments)
        )


class SingleRefProfile(object):
    """
    The SingleRefProfile class encapsulates the alignment profile of a single reference sequence.

    This class stores the length of the reference sequence, the number of replicates, the length of sRNA, and a list of all sRNA alignments against this reference.

    Attributes:
    ref_len : int
        The length of the reference sequence.
    replicates : int
        The number of replicates.
    srna_len : int
        The length of sRNA.
    all_alignments : list of SingleAlignment objects
        A list of SingleAlignment objects representing all sRNA alignments against the reference sequence.

    Methods:
    __str__():
        Returns a string representation of the SingleRefProfile object, suitable for output.
    __repr__():
        Returns the same as __str__(), a string representation of the object.
    __eq__(other: SingleRefProfile) -> bool:
        Determines if two SingleRefProfile objects are identical, including their list of alignments.
        `other` is the other SingleRefProfile object to compare with.
    """

    def __init__(self):
        self.ref_len = 0
        self.replicates = 0
        self.srna_len = 0
        # Columnar storage - one entry per alignment. SingleAlignment objects are
        # only materialized on demand (all_alignments property), so parsing a
        # chromosome-scale CSV does not allocate millions of small objects.
        self.srnas = []
        self.positions = None  # int64 array
        self.strands = None  # bool array, True for "+"
        self.times = None  # int64 array
        self.counts = None  # (n, replicates) float array
        self._alignments_cache = None
        self._alignments_override = None

    @property
    def all_alignments(self):
        if self._alignments_override is not None:
            return self._alignments_override
        if self.positions is None:
            return []
        if self._alignments_cache is None:
            self._alignments_cache = [
                SingleAlignment(DNA(s), int(p), "+" if f else "-", int(t), c)
                for s, p, f, t, c in zip(
                    self.srnas, self.positions, self.strands, self.times, self.counts
                )
            ]
        return self._alignments_cache

    @all_alignments.setter
    def all_alignments(self, alignments):
        try:
            self._set_columns_from(alignments)
            self._alignments_override = None
        except (AttributeError, TypeError):
            # not SingleAlignment-like objects; keep the raw list (used by tests)
            self._alignments_override = list(alignments)

    def _set_columns_from(self, alignments):
        alignments = list(alignments)
        if not alignments:
            self.srnas = []
            self.positions = None
            self.strands = None
            self.times = None
            self.counts = None
            self._alignments_cache = None
            return
        positions = np.array([sa.position for sa in alignments], dtype=np.int64)
        strands = np.array([sa.strand == "+" for sa in alignments])
        times = np.array([sa.times_aligned for sa in alignments], dtype=np.int64)
        counts = np.stack(
            [np.asarray(sa.indv_alignments, dtype=float) for sa in alignments]
        )
        self.srnas = [str(sa.srna) for sa in alignments]
        self.positions = positions
        self.strands = strands
        self.times = times
        self.counts = counts
        self._alignments_cache = alignments

    def __str__(self):
        return "{0}\t{1}".format(self.ref_len, self.all_alignments)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.ref_len == other.ref_len
            and self.all_alignments == other.all_alignments
            and self.replicates == other.replicates
            and self.srna_len == other.srna_len
        )


class RefProfiles(object):
    """
    The RefProfiles class represents the alignment profiles for all reference sequences in an alignment file of a set read or sRNA length.

    This class stores the length of the small RNA read (sRNA), the number of replicates, and a dictionary containing SingleRefProfile
    objects for each reference sequence.

    Attributes:
    srna_len : int
        The length of the sRNA.
    replicates : int
        The number of replicate samples.
    single_ref_profiles : dict
        A dictionary where keys are reference sequence headers and values are SingleRefProfile objects representing the alignment profile for each reference sequence.

    Methods:
    __str__():
        Returns a string representation of the RefProfiles object, suitable for output.
    __eq__():
        Determines if two RefProfiles objects are identical.
    load_single_ref_profiles(in_file):
        Loads alignment profiles from a given scram2 profile file. Each line of the file represents a SingleAlignment of a sRNA read.
    """

    def __init__(self):
        self.srna_len = 0
        self.replicates = 0
        self.single_ref_profiles = {}
        self.subset = False

    def __str__(self):
        return "{0}\t{1}\t{2}".format(
            self.srna_len, self.replicates, self.single_ref_profiles
        )

    def __eq__(self, other):
        return (
            self.replicates == other.replicates
            and self.single_ref_profiles == other.single_ref_profiles
            and self.srna_len == other.srna_len
        )

    def load_single_ref_profiles(
        self, in_file, header=None, start=None, end=None, padding=0
    ):
        """Loads a scram2 profile file for a single sRNA length

        Args:
            in_file (string): scram2 profile file path
            header (string): specific header to be parsed
            start (int): start position of the range
            end (int): end position of the range
            padding (int): number of base pairs to include on either end of the range
        """
        subset = start is not None and end is not None
        self.subset = subset

        df = pd.read_csv(in_file)
        columns = [str(c) for c in df.columns]
        if "Std. Err" in columns:
            raise ValueError(
                f"{in_file} is in the aligner's mean/standard-error format; "
                "scram2Plot requires the per-replicate format "
                "(Header,len,sRNA,Position,Strand,Times aligned,<file columns>)"
            )
        if len(columns) < 7 or columns[0] != "Header":
            raise ValueError(
                f"{in_file} is not a scram profile CSV; expected columns "
                "Header,len,sRNA,Position,Strand,Times aligned,<file columns>"
            )
        replicates = len(columns) - 6

        headers = df.iloc[:, 0].astype(str).str.split().str[0].to_numpy()
        ref_lens = df.iloc[:, 1].to_numpy(dtype=np.int64)
        srnas = df.iloc[:, 2].astype(str).to_numpy()
        positions = df.iloc[:, 3].to_numpy(dtype=np.int64)
        strands = df.iloc[:, 4].to_numpy() == "+"
        times = df.iloc[:, 5].to_numpy(dtype=np.int64)
        counts = df.iloc[:, 6:].to_numpy(dtype=float)
        if np.isnan(counts).any():
            raise ValueError(
                f"{in_file} has missing or non-numeric read-count values"
            )

        # Keep only rows for the requested header and window. A read overlapping
        # the left edge of the window still covers window positions, so the left
        # cut-off is widened by each read's length; DataForPlot clips the
        # overhang after coverage conversion.
        keep = np.ones(len(df), dtype=bool)
        if header is not None:
            keep &= headers == header
        if subset:
            srna_lens = df.iloc[:, 2].astype(str).str.len().to_numpy(dtype=np.int64)
            keep &= positions + srna_lens - 1 >= start - padding
            keep &= positions <= end + padding
        if not keep.all():
            headers = headers[keep]
            ref_lens = ref_lens[keep]
            srnas = srnas[keep]
            positions = positions[keep]
            strands = strands[keep]
            times = times[keep]
            counts = counts[keep]
        if subset:
            # positions become relative to the subset start (possibly negative
            # for reads overlapping the left edge)
            positions = positions - start + padding

        if len(headers) == 0:
            raise ValueError(
                f"{in_file} contains no alignments"
                + (f" for header {header}" if header is not None else "")
            )

        group_codes, unique_headers = pd.factorize(headers)
        for code, row_header in enumerate(unique_headers):
            mask = group_codes == code
            profile = SingleRefProfile()
            profile.ref_len = (
                end - start + 2 * padding if subset else int(ref_lens[mask][0])
            )
            profile.replicates = replicates
            profile.srnas = srnas[mask].tolist()
            profile.positions = positions[mask]
            profile.strands = strands[mask]
            profile.times = times[mask]
            profile.counts = counts[mask]
            profile.srna_len = len(profile.srnas[0])
            self.single_ref_profiles[row_header] = profile
        self.srna_len = self.single_ref_profiles[unique_headers[0]].srna_len
        self.replicates = replicates


class SectStruct(object):
    """
    Class for handling secondary structure data.

    Attributes:
    - ref_file: The path to the reference FASTA file containing sequences.
    - header: The identifier of the sequence to use from the FASTA file.
    - start: Optional start position of the sequence. Default is None.
    - end: Optional end position of the sequence. Default is None.
    - sequence: The subsequence from start to end from the FASTA file.
    - structure: The secondary structure of the sequence.
    - min_free_energy: minimum free energy of the structure
    - pairs: A list of base-pair positions in the structure.
    - length: The length of the sequence.
    """

    def __init__(self, ref_file, header, start=None, end=None):
        """
        Initializes a new SectStruct object.

        Parameters:
        - ref_file: The path to the reference FASTA file.
        - header: Identifier for the sequence of interest in the FASTA file.
        - start: Optional start position (0-based index) for the sequence. Default is None.
        - end: Optional end position (0-based index) for the sequence. Default is None.
        """
        self.ref_file = ref_file
        self.header = header
        self.start = start
        self.end = end
        self.sequence = None
        self.structure = None
        self.min_free_energy = 0.0
        self.pairs = []

        # Initialize the sequence, structure, and base pairs.
        self._get_sequence()
        self.length = len(self.sequence)
        self._get_structure()
        self._get_pairs()

    def _get_sequence(self):
        """
        Fetches and stores the sequence based on the given header and optionally start and end positions.

        Raises:
        - FileNotFoundError if the reference file is not found.
        - ValueError if the sequence is not found in the file or if start/end positions are invalid.
        """
        from Bio import SeqIO  # deferred: only needed for secondary-structure plots

        # Check if the reference file exists
        if not os.path.isfile(self.ref_file):
            raise FileNotFoundError(f"Reference file {self.ref_file} not found.")

        sequence_found = False
        for seq_record in SeqIO.parse(self.ref_file, "fasta"):
            if seq_record.id == self.header:
                sequence_found = True

                if self.start is not None and self.end is not None:
                    # Check bounds and order of start and end positions
                    if self.start < 0 or self.end > len(seq_record.seq):
                        raise ValueError(
                            f"Start and end positions must be within the bounds of the sequence (0, {len(seq_record.seq)})"
                        )
                    if self.end <= self.start:
                        raise ValueError(f"End position must be after start position")

                    self.sequence = str(seq_record.seq[self.start : self.end])
                else:
                    self.sequence = str(seq_record.seq)

                break

        if not sequence_found:
            raise ValueError(
                f"Sequence with header {self.header} not found in reference file."
            )

    def _get_structure(self):
        """
        Computes the secondary structure and minimum free energy of the stored sequence using RNA.fold and stores it.
        """
        from ViennaRNA import RNA  # deferred: only needed for secondary-structure plots

        self.structure, self.min_free_energy = RNA.fold(self.sequence)

    def _get_pairs(self):
        """
        Calculates and stores the base-pairing positions from the secondary structure.
        """
        stack = []
        for i, char in enumerate(self.structure):
            if char == "(":
                stack.append(i + 1)  # index starts from 1
            elif char == ")":
                start = stack.pop()
                end = i + 1  # index starts from 1
                self.pairs.append((start, end))


class DataForPlot(object):
    """
    Class for managing data relevant to plotting biological sequence data.

    Attributes:
    - ref_profiles: A reference profiles object containing multiple sequence alignments.
    - header: Header string to identify a specific sequence.
    - ref_len: Length of the reference sequence specified by header.
    - srna_len: Length of sRNA.
    - replicates: Number of replicate data sets.
    - fwd: Array containing forward strand alignment data.
    - rvs: Array containing reverse strand alignment data.
    - x_axis: List of x-axis values for plotting.
    - y_flat: Flattened list of y-axis values for plotting.
    """

    def __init__(self, ref_profiles, header, abund=False):
        """
        Initialize a new DataForPlot object.

        Args:
        - ref_profiles: A reference profiles object containing sequence data.
        - header: Header string to specify which sequence to use.
        - abund: Whether to extract data using abundance (max) or coverage (sum) logic.
        """
        self.header = header
        profile = ref_profiles.single_ref_profiles[header]
        self.ref_len = profile.ref_len
        self.srna_len = ref_profiles.srna_len
        self.replicates = ref_profiles.replicates
        # For subset windows, reads overlapping the left edge have negative
        # relative positions; extra leading rows hold them until coverage has
        # been smeared, then trim_margin() removes the overhang.
        self.margin = self.srna_len - 1 if ref_profiles.subset else 0
        rows = self.margin + self.ref_len + 1
        self.fwd = np.zeros((rows, self.replicates), dtype=float)
        self.rvs = np.zeros((rows, self.replicates), dtype=float)
        self.x_axis = np.arange(self.ref_len + 1)
        self.y_flat = []
        self._extract_from_ref_profiles(profile, abund)

    def _extract_from_ref_profiles(self, profile, abund):
        """Extracts alignment data for a specific sequence from the profile's columns."""
        if profile.positions is None or len(profile.positions) == 0:
            return
        positions = profile.positions + self.margin
        counts = profile.counts
        for arr, mask in ((self.fwd, profile.strands), (self.rvs, ~profile.strands)):
            p = positions[mask]
            c = counts[mask]
            if abund:
                # highest-mean read wins each position: assign in ascending mean
                # order so the max lands last
                order = np.argsort(c.mean(axis=1), kind="stable")
                arr[p[order]] = c[order]
            else:
                np.add.at(arr, p, c)

    def trim_margin(self):
        """Removes the left-overhang rows used for reads spanning the left edge of a
        subset window (a no-op for whole-reference data)."""
        if self.margin:
            self.fwd = self.fwd[self.margin :]
            self.rvs = self.rvs[self.margin :]
            self.margin = 0

    def convert_to_coverage(self, abund=True):
        """Converts raw alignment data to coverage data."""
        self.fwd = self._coverage_per_strand(self.fwd, abund)
        self.rvs = self._coverage_per_strand(self.rvs, abund)
        self.trim_margin()

    def _coverage_per_strand(self, old_array, abund):
        """
        Convert to coverage per strand: each read's counts are spread across the
        srna_len positions its alignment covers.

        Args:
        - old_array: Existing array of alignment data (reads at their 5'-most position).
        - abund: If True, each position shows the counts of the highest-mean read
          covering it (winner takes all); otherwise counts of overlapping reads sum.

        Returns:
        - new_arr: A new array containing coverage data.
        """
        n, reps = old_array.shape
        window_starts = n - self.srna_len + 1  # valid alignment start rows
        if window_starts <= 0:
            return np.zeros_like(old_array)

        if abund:
            means = old_array[:window_starts].mean(axis=1)
            pad = np.full(self.srna_len - 1, -np.inf)
            padded = np.concatenate((pad, means, pad))
            # windows[j] holds the means of every read window covering position j
            windows = np.lib.stride_tricks.sliding_window_view(padded, self.srna_len)[:n]
            best = windows.max(axis=1)
            # np.argmax takes the first max, matching the original first-window-wins
            # tie handling
            src = windows.argmax(axis=1) + np.arange(n) - (self.srna_len - 1)
            return np.where(
                (best > 0)[:, None],
                old_array[np.clip(src, 0, window_starts - 1)],
                0.0,
            )

        diff = np.zeros((n + self.srna_len, reps), dtype=old_array.dtype)
        diff[:window_starts] = old_array[:window_starts]
        diff[self.srna_len : window_starts + self.srna_len] -= old_array[:window_starts]
        return np.cumsum(diff, axis=0)[:n]

    def convert_to_error_bounds(self):
        """Calculates and stores error bounds for the data."""
        self.fwd = self._error_bounds(self.fwd)
        self.rvs = self._error_bounds(self.rvs)

    def _error_bounds(self, old_array):
        """
        Calculate mean +/- standard error bounds for an array.

        Uses the sample standard deviation (ddof=1), as required for the standard
        error of a mean across replicates. With a single replicate the error is 0.

        Args:
        - old_array: Existing array of data.

        Returns:
        - new_arr: A (rows, 2) array containing lower and upper error bounds.
        """
        mean = old_array.mean(axis=1)
        reps = old_array.shape[1]
        if reps > 1:
            se = old_array.std(axis=1, ddof=1) / np.sqrt(reps)
        else:
            se = np.zeros_like(mean)
        return np.column_stack((mean - se, mean + se))

    def flatten(self, d=1):
        """
        Flattens data arrays for both strands into self.y_flat: all forward-strand
        columns first, then all (negated) reverse-strand columns, so y_flat[i] and
        y_flat[i + n_columns] are the same column's two strands. Replaces any
        previous contents of y_flat.

        Args:
        - d: Degree of smoothing. Default is 1 (no smoothing).
        """
        self.y_flat = []
        for i in range(self.fwd.shape[1]):
            self.y_flat.append(self._maybe_smooth(self.fwd[:, i], d))
        for i in range(self.rvs.shape[1]):
            self.y_flat.append(self._maybe_smooth(-self.rvs[:, i], d))

    @classmethod
    def _maybe_smooth(cls, column, d):
        return cls._smoothTriangle(column, d) if d > 1 else column

    @staticmethod
    def _smoothTriangle(data, degree):
        """
        Smooths data using a centred triangular window; partial windows at the
        edges are renormalized rather than padded with repeated values.

        Args:
        - data: List or array of data points.
        - degree: Degree of smoothing.

        Returns:
        - smoothed: Smoothed data as a numpy array.
        """
        data = np.asarray(data, dtype=float)
        if degree < 1 or len(data) == 0:
            return data
        kernel = np.concatenate(
            (np.arange(degree + 1), np.arange(degree)[::-1])
        ).astype(float)  # up then down, centred, length 2*degree+1
        kernel /= kernel.sum()
        smoothed = np.convolve(data, kernel, mode="same")
        edge_weight = np.convolve(np.ones_like(data), kernel, mode="same")
        return smoothed / edge_weight

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.header,
            self.ref_len,
            self.srna_len,
            self.fwd,
            self.rvs,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.header == other.header
            and self.ref_len == other.ref_len
            and self.srna_len == other.srna_len
            and np.array_equal(self.fwd, other.fwd)
            and np.array_equal(self.rvs, other.rvs)
        )


class AlignmentPlot:
    """
    A class for generating alignment plots of biological sequence data.

    Attributes:
        align_prefix (str): Prefix for the alignment files.
        align_lens (list): List of alignment lengths.
        header (str or list): Header for the sequence or list of headers.
        start (int, optional): Start position for the sequence. Defaults to None.
        end (int, optional): End position for the sequence. Defaults to None.
        smoothing_window (int, optional): Window size for smoothing. Defaults to 1.
        cov (bool, optional): Whether to convert to coverage. Defaults to True.
        abund (bool, optional): Whether to include abundance data. Defaults to True.
        se (bool, optional): Whether to include standard error bounds. Defaults to True.
        save (bool, optional): Whether to save the plot. Defaults to True.
        ylim_set (float, optional): Set y-axis limits. Defaults to 0.
        sec_structure (bool, optional): Whether to include secondary structure. Defaults to False.
        ref_file (str, optional): Reference file for secondary structure. Defaults to None.
        show_seq (bool, optional): Whether to show the sequence. Defaults to False.
    """

    def __init__(
        self,
        align_prefix,
        align_lens,
        header,
        start=None,
        end=None,
        smoothing_window=1,
        cov=True,
        abund=True,
        se=True,
        save=True,
        ylim_set=0,
        sec_structure=False,
        ref_file=None,
        show_seq=False,
    ):
        self.align_prefix = align_prefix
        self.align_lens = align_lens
        self.header = header
        self.start = start
        self.end = end
        self.smoothing_window = smoothing_window
        self.cov = cov
        self.abund = abund
        self.se = se
        self.save = save
        self.ylim_set = ylim_set
        self.sec_structure = sec_structure
        self.ref_file = ref_file
        self.show_seq = show_seq
        self.ss = None
        self.all_handles = []
        self.all_labels = []

    def set_up_plot(self):
        """
        Sets up the current matplotlib figure for displaying biological sequence data.
        """
        # Set the x-axis label
        plt.xlabel("Position")

        # If a secondary structure object is provided, perform additional setup
        # (sec_struct_setup applies its own equal aspect so the arcs stay circular)
        if self.ss is not None:
            self.sec_struct_setup()

        # Hide y-axis ticks
        plt.gca().set_yticks([])

        # Get the current axis and set the y-label
        ax = plt.gca()
        ax.set_ylabel("Abundance")

        # Move the y-label to the left
        ax.yaxis.set_label_coords(-0.04, 0.5)

    def sec_struct_setup(self):
        """
        Adds secondary structure information to an existing matplotlib plot.
        """
        max_y = 0  # Variable to store maximum y value

        for x1, x2 in self.ss.pairs:
            # Define the two points
            p1 = np.array([x1, 0])
            p2 = np.array([x2, 0])

            # Calculate the center of the circle
            c = (p1 + p2) / 2

            # Calculate the radius of the circle
            r = np.linalg.norm(p1 - p2) / 2
            # Update max_y if the current radius is greater
            max_y = max(max_y, r)

            # Define the angles for the arc
            t = np.linspace(0, np.pi, 100)

            # Parametric equations for the circle
            x = c[0] + r * np.cos(t)
            y = c[1] + r * np.sin(t)

            # Plot the arc
            plt.plot(x, y, c="blue", linewidth=0.1, alpha=0.8, antialiased=True)

        # Ensure the circle isn't distorted
        plt.axis("equal")

        # Set the x-axis limits based on sequence length
        plt.xlim(0, self.ss.length + 1)

        # Set symmetrical y-limits around 0
        plt.ylim(-max_y, max_y)

    def generate_plot(self):
        """
        Generate the alignment plot based on the provided configurations.
        """
        if self.sec_structure:
            self.ss = SectStruct(self.ref_file, self.header, self.start, self.end)
            if self.show_seq:
                print("Selected sequence:")
                print(self.ss.sequence)

        prepared = self._prepare_data()
        if not prepared:
            raise RuntimeError(
                f"No alignments found for header {self.header} at lengths "
                f"{self.align_lens} with prefix {self.align_prefix}"
            )

        fig = plt.figure(figsize=(12, 12))
        try:
            self.set_up_plot()
            if self.ylim_set == 0:
                self.ylim_set = self._max_abundance(prepared) * 1.1
            for spd in prepared.values():
                self.single_plot(spd, self.header)
            if self.all_handles:
                plt.legend(self.all_handles, self.all_labels, loc="upper right")
            if self.save:
                self.save_plot()
            plt.show()
        finally:
            plt.close(fig)

    def _prepare_data(self):
        """
        Loads and transforms the alignment data exactly once per sRNA length.

        Returns:
            dict: {srna_len: DataForPlot} for every length with alignments.
        """
        prepared = {}
        for srna_len in self.align_lens:
            file_path = f"{self.align_prefix}_{srna_len}.csv"
            if not os.path.isfile(file_path):
                print(f"File {file_path} not found. Skipping.")
                continue
            rp = RefProfiles()
            try:
                rp.load_single_ref_profiles(
                    file_path, header=self.header, start=self.start, end=self.end
                )
            except ValueError as e:
                print(f"{e}. Skipping.")
                continue
            spd = DataForPlot(rp, self.header, abund=self.abund)
            if self.cov:
                spd.convert_to_coverage(abund=self.abund)
            else:
                spd.trim_margin()
            if self.se:
                spd.convert_to_error_bounds()
            spd.flatten(self.smoothing_window)
            prepared[srna_len] = spd
        return prepared

    @staticmethod
    def _max_abundance(prepared):
        """
        The largest absolute plotted value across all prepared lengths, both strands.
        """
        y_max = 0.0
        for spd in prepared.values():
            n_cols = spd.fwd.shape[1]
            for i in range(n_cols):
                y_max = max(
                    y_max, spd.y_flat[i].max(), abs(spd.y_flat[i + n_cols].min())
                )
        return y_max

    def single_plot(self, spd, header):
        """
        Plots the prepared data for a single sRNA length.

        Args:
            spd (DataForPlot): Prepared (transformed and flattened) plot data.
            header (str): The header of the sequence.
        """
        ax1 = plt.gca()  # primary axis
        ax2 = ax1.twinx()  # secondary axis

        self.configure_plot_axes(ax1, ax2)
        ax2.yaxis.tick_left()

        ax2.plot(spd.x_axis, np.zeros(len(spd.x_axis)), color="grey", linewidth=0.5)
        self.set_plot_title(ax2, header)

        color = color_for_length(spd.srna_len)
        if self.se:
            self.plot_with_error_bounds(ax2, spd, color)
        else:
            self.plot_reps_individually(ax2, spd, color)

        if self.start is not None:
            x_ticks = ax2.get_xticks()
            ax2.set_xticks(x_ticks, [int(x + self.start) for x in x_ticks])

        ax2.yaxis.set_label_coords(-0.1, 0.5)
        handles, labels = ax2.get_legend_handles_labels()
        self.all_handles.extend(handles)
        self.all_labels.extend(labels)
        return handles, labels

    def configure_plot_axes(self, ax1, ax2):
        """
        Configures the primary and secondary axes of the plot.

        Args:
            ax1: The primary axis.
            ax2: The secondary axis.
        """
        ax1.yaxis.set_ticklabels([])
        ax1.yaxis.set_ticks([])
        ax2.set_ylim(-self.ylim_set, self.ylim_set)
        ax2.yaxis.tick_left()

    def set_plot_title(self, ax2, header):
        """
        Sets the title of the plot.

        Args:
            ax2: The secondary axis.
            header (str): The header of the sequence.
        """
        header = header.split()[0]
        title = f"Profile for {header}"
        if self.start is not None and self.end is not None:
            title += f" from position {self.start} to {self.end}"
        ax2.set_title(title)

    def plot_with_error_bounds(self, ax2, spd, color):
        """
        Plots the mean +/- standard error ribbon for each strand.

        Args:
            ax2: The secondary axis.
            spd (DataForPlot): Prepared plot data; y_flat holds
                [fwd_lower, fwd_upper, rvs_lower, rvs_upper].
            color: Plot color for this sRNA length.
        """
        ax2.fill_between(
            spd.x_axis,
            spd.y_flat[0],
            spd.y_flat[1],
            color=color,
            alpha=0.4,
            label=str(spd.srna_len) + " nt",
        )
        ax2.fill_between(
            spd.x_axis,
            spd.y_flat[2],
            spd.y_flat[3],
            color=color,
            alpha=0.4,
        )

    def plot_reps_individually(self, ax2, spd, color):
        """
        Plots each replicate's forward and reverse strand traces.

        Args:
            ax2: The secondary axis.
            spd (DataForPlot): Prepared plot data; y_flat holds all forward
                replicate columns followed by all (negated) reverse columns.
            color: Plot color for this sRNA length.
        """
        for i in range(spd.replicates):
            fwd = spd.y_flat[i]
            rvs = spd.y_flat[i + spd.replicates]
            label = str(spd.srna_len) + " nt" if i == 0 else None
            ax2.plot(spd.x_axis, fwd, color=color, alpha=0.8, label=label)
            ax2.plot(spd.x_axis, rvs, color=color, alpha=0.8)
            ax2.fill_between(spd.x_axis, fwd, rvs, color=color, alpha=0.05)

    def save_plot(self):
        if isinstance(self.header, list):
            save_file = self.align_prefix + "_" + "_".join(self.header) + ".svg"
        else:
            save_file = self.align_prefix + "_" + self.header + ".svg"

        if self.start is not None and self.end is not None:
            save_file = save_file.replace(".svg", f"_{self.start}_{self.end}.svg")

        plt.savefig(save_file)


def comma_separated_ints(value):
    """
    Parses a comma-separated string of integers into a list.
    """
    try:
        return [int(v) for v in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Invalid comma-separated integers: '{}'".format(value)
        )


def comma_separated_strings(value):
    return value.split(",")


def main():
    """
    Entry point of the script. Parses the command-line arguments and triggers plotting.

    Notes:
        - Calls `align_plot` function with the parsed arguments.
    """

    # Argument parser
    parser = argparse.ArgumentParser(description="Plot abundance profiles")

    # Required arguments
    parser.add_argument("align_prefix", help="Prefix of alignment files")
    parser.add_argument(
        "align_lens",
        help="Comma-separated list of siRNA lengths to plot",
        type=comma_separated_ints,
    )
    parser.add_argument("header", help="Reference header to plot", type=str)

    # Optional arguments
    parser.add_argument(
        "-s", "--smoothing_window", help="Smoothing window", type=int, default=1
    )
    parser.add_argument("-c", "--coverage", help="Plot coverage", action="store_true")
    parser.add_argument("-a", "--abundance", help="Plot abundance", action="store_true")
    parser.add_argument("-e", "--error", help="Plot error", action="store_true")
    parser.add_argument(
        "-y", "--ylim", help="Set y-axis limit", type=int, default=0
    )
    parser.add_argument("-n", "--no_save", help="Do not save plot", action="store_true")
    parser.add_argument(
        "-start", "--start", help="Start position for the subset", type=int
    )
    parser.add_argument("-end", "--end", help="End position for the subset", type=int)
    parser.add_argument(
        "--sec_structure",
        help="Plot secondary structure of reference sequence",
        action="store_true",
    )
    parser.add_argument(
        "--ref_file",
        help="Alignment reference file (FASTA) for use in generating secondary structure",
        type=str,
        default=None,
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Error checking: Required ref_file for sec_structure
    if args.sec_structure and args.ref_file is None:
        parser.error("--ref_file is required when --sec_structure is set")

    # If abundance is specified, enable coverage as well
    if args.abundance:
        args.coverage = True

    # Create an instance of AlignmentPlot
    plotter = AlignmentPlot(
        align_prefix=args.align_prefix,
        align_lens=args.align_lens,
        header=args.header,
        start=args.start,
        end=args.end,
        smoothing_window=args.smoothing_window,
        cov=args.coverage,
        abund=args.abundance,
        se=args.error,
        save=not args.no_save,
        ylim_set=args.ylim,
        sec_structure=args.sec_structure,
        ref_file=args.ref_file,
    )

    # Generate the plot
    try:
        plotter.generate_plot()
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


# Note: Add the following lines to trigger the main() function when the script is executed:
if __name__ == "__main__":
    main()
