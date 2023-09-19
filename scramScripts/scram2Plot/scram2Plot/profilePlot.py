import numpy as np
import matplotlib.pyplot as plt
import os.path
import csv
import argparse
from Bio import SeqIO
import RNA




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
        self.all_alignments = []

    def __str__(self):
        return "{0}\t{1}".format(self.ref_len, self.all_alignments)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.ref_len == other.ref_len
            and self.all_alignments == other.all_alignments
            and self.replicates == other.replicates
            and self.snra_len == other.srna_len
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
        try:
            with open(in_file, "r") as in_handle:
                reader = csv.reader(in_handle, delimiter=",")
                for row in reader:
                    if row[0] == "Header":
                        continue
                    (
                        row_header,
                        ref_len,
                        srna,
                        position,
                        strand,
                        times_aligned,
                        *indv_alignments,
                    ) = row
                    row_header = row_header.split()[0]
                    # Skip alignments that do not match the input header or are outside the input positions
                    if header is not None and row_header != header:
                        continue
                    if (
                        start is not None
                        and end is not None
                        and (
                            int(position) < start - padding
                            or int(position) > end + padding
                        )
                    ):
                        continue

                    # Construct a SingleAlignment object

                    srna = DNA(srna)
                    times_aligned = int(times_aligned)
                    ref_len = int(ref_len)
                    position = int(position)
                    if start is not None:
                        position = (
                            position - start + padding
                        )  # New position relative to the start of the subset
                        ref_len = end - start + 2 * padding  # New reference length
                    indv_alignments = np.array([float(x) for x in indv_alignments])
                    sa = SingleAlignment(
                        srna, position, strand, times_aligned, indv_alignments
                    )

                    # Add the SingleAlignment to the corresponding SingleRefProfile
                    if row_header not in self.single_ref_profiles:
                        self.single_ref_profiles[row_header] = SingleRefProfile()
                        self.single_ref_profiles[row_header].ref_len = ref_len
                        self.single_ref_profiles[row_header].replicates = len(
                            indv_alignments
                        )
                        self.single_ref_profiles[row_header].srna_len = len(srna)
                        self.srna_len = len(srna)
                    self.single_ref_profiles[row_header].all_alignments.append(sa)

                # Set the number of replicates based on the last read alignment
                if "sa" in locals():
                    self.replicates = len(sa.indv_alignments)
                else:
                    raise ValueError("The input file does not contain any alignments")
        except FileNotFoundError:
            print(f"The input file {in_file} does not exist. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing the input file {in_file}: {e}")


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

    def __init__(self, ref_profiles, header):
        """
        Initialize a new DataForPlot object.

        Args:
        - ref_profiles: A reference profiles object containing sequence data.
        - header: Header string to specify which sequence to use.
        """
        self.ref_profiles = ref_profiles
        self.header = header
        self.ref_len = ref_profiles.single_ref_profiles[header].ref_len
        self.srna_len = ref_profiles.srna_len
        self.replicates = ref_profiles.replicates
        self.fwd = np.zeros((self.ref_len + 1, self.replicates), dtype=float)
        self.rvs = np.zeros((self.ref_len + 1, self.replicates), dtype=float)
        self.x_axis = list(range(len(self.fwd)))
        self.y_flat = []
        self._extract_from_ref_profiles()

    def _extract_from_ref_profiles(self):
        """Extracts alignment data for a specific sequence from the reference profiles object."""
        for sa in self.ref_profiles.single_ref_profiles[self.header].all_alignments:
            if sa.strand == "+":
                self.fwd[sa.position] = sa.indv_alignments
            else:
                self.rvs[sa.position] = sa.indv_alignments

    def convert_to_coverage(self, abund=True):
        """Converts raw alignment data to coverage data."""
        self.fwd = self._coverage_per_strand(self.fwd, abund)
        self.rvs = self._coverage_per_strand(self.rvs, abund)

    def _coverage_per_strand(self, old_array, abund):
        """
        Convert to coverage per strand.

        Args:
        - old_array: Existing array of alignment data.
        - abund: If True, abundance-based calculation will be performed.

        Returns:
        - new_arr: A new array containing coverage data.
        """
        new_arr = np.zeros((old_array.shape), dtype=float)
        for i in range(len(new_arr) - self.srna_len + 1):
            for j in range(len(new_arr[i])):
                if not abund:
                    new_arr[i : i + self.srna_len, j] += old_array[i, j]
                else:
                    for k in range(self.srna_len):
                        if old_array[i, j] > new_arr[i + k, j]:
                            new_arr[i + k, j] = old_array[i, j]
        return new_arr

    def convert_to_error_bounds(self):
        """Calculates and stores error bounds for the data."""
        self.fwd = self._error_bounds(self.fwd)
        self.rvs = self._error_bounds(self.rvs)

    def _error_bounds(self, old_array):
        """
        Calculate error bounds for an array.

        Args:
        - old_array: Existing array of data.

        Returns:
        - new_arr: A new array containing lower and upper error bounds.
        """
        new_arr = np.zeros((self.ref_len + 1, 2), dtype=float)
        for i in range(len(new_arr)):
            new_arr[i, 0] = np.mean(old_array[i, :]) - (
                np.std(old_array[i, :]) / np.sqrt(self.replicates)
            )
            new_arr[i, 1] = np.mean(old_array[i, :]) + (
                np.std(old_array[i, :]) / np.sqrt(self.replicates)
            )
        return new_arr

    def flatten(self, d=1):
        """
        Flattens data arrays for both strands and appends them to self.y_flat.

        Args:
        - d: Degree of smoothing. Default is 1.
        """
        if d == 1:
            for i in range(self.fwd.shape[1]):
                self.y_flat.append(self.fwd[:, [i]].flatten())
                self.y_flat.append(-self.rvs[:, [i]].flatten())
        else:
            for i in range(self.fwd.shape[1]):
                self.y_flat.append(self._smoothTriangle(self.fwd[:, [i]].flatten(), d))
                self.y_flat.append(self._smoothTriangle(-self.rvs[:, [i]].flatten(), d))

    @staticmethod
    def _smoothTriangle(data, degree):
        """
        Smooths data using a triangular smoothing function.

        Args:
        - data: List or array of data points.
        - degree: Degree of smoothing.

        Returns:
        - smoothed: Smoothed data.
        """
        triangle = np.concatenate(
            (np.arange(degree + 1), np.arange(degree)[::-1])
        )  # up then down
        smoothed = []

        for i in range(degree, len(data) - degree * 2):
            point = data[i : i + len(triangle)] * triangle
            smoothed.append(np.sum(point) / np.sum(triangle))
        # Handle boundaries
        smoothed = [smoothed[0]] * int(
            degree + degree / 2
        ) + smoothed  # TODO: this can be better
        while len(smoothed) < len(data):
            smoothed.append(smoothed[-1])
        return smoothed

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


def set_up_plot(ss=None):
    """
    Sets up a matplotlib plot for displaying biological sequence data.

    Args:
        ss (SectStruct, optional): A SectStruct object containing information about the secondary structure.
                                   If provided, additional setup related to the secondary structure will be performed.

    Notes:
        - Creates a 12x12-inch plot with a DPI of 600.
        - X-axis label is set to "Position".
        - Y-axis label is set to "Abundance", and it is placed slightly to the left of the axis.
        - Y-axis ticks are hidden.
        - If a SectStruct object is provided, additional setup is performed via the `sec_struct_setup` function (not defined here).
        - TODO: Verify if setting the axis to "equal" is necessary for your use case.
    """

    # Create a new figure with specified size and DPI
    plt.figure(figsize=(12, 12), dpi=600)

    # Set the x-axis label
    plt.xlabel("Position")

    # Uncomment if you decide to keep the y-axis label as "Abundance"
    # plt.ylabel("Abundance")

    # Make the x and y dimensions equal; TODO: Check if this is needed
    plt.axis("equal")

    # If a secondary structure object is provided, perform additional setup
    if ss is not None:
        sec_struct_setup(ss)  # Note: `sec_struct_setup` needs to be defined elsewhere

    # Hide y-axis ticks
    plt.gca().set_yticks([])

    # Get the current axis and set the y-label
    ax = plt.gca()
    ax.set_ylabel("Abundance")

    # Move the y-label to the left
    ax.yaxis.set_label_coords(-0.04, 0.5)


def sec_struct_setup(ss):
    """
    Adds secondary structure information to an existing matplotlib plot.

    Args:
        ss (SectStruct): A SectStruct object that contains information about
                         the secondary structure of the sequence.
                         Expected to have attributes 'pairs' and 'length'.

    Notes:
        - Plots arcs between paired nucleotide positions in the sequence.
        - The arcs are plotted in blue with a linewidth of 0.1 and alpha of 0.8.
        - The x-axis limits are set based on the sequence length.
        - The y-axis limits are set symmetrically around 0, based on the maximum radius of the arcs.
        - The aspect ratio is set to 'equal' to avoid distortion.
    """

    max_y = 0  # Variable to store maximum y value

    for x1, x2 in ss.pairs:
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
    plt.xlim(0, ss.length + 1)

    # Set symmetrical y-limits around 0
    plt.ylim(-max_y, max_y)



def single_plot(
    ref_profiles,
    header,
    start=None,
    end=None,
    smoothing_window=1,
    cov=True,
    abund=True,
    se=True,
    ylim_set=0,
    padding=0,
):
    """
    Plots the data for a single sequence based on its reference profiles.

    Args:
        ref_profiles (RefProfiles): The reference profiles object containing alignment and coverage data.
        header (str): The header of the reference sequence.
        start (int, optional): The start position in the sequence for the plot.
        end (int, optional): The end position in the sequence for the plot.
        smoothing_window (int, optional): The size of the smoothing window. Default is 1 (no smoothing).
        cov (bool, optional): Whether to include coverage in the plot. Default is True.
        abund (bool, optional): Whether to include abundance in the plot. Default is True.
        se (bool, optional): Whether to include standard error bounds in the plot. Default is True.
        ylim_set (int, optional): The y-limit for the plot. Default is 0 (automatically determined).
        padding (int, optional): Padding to adjust the x-axis ticks. Default is 0.

    Returns:
        tuple: A tuple containing handles and labels for the plot legend, if successful.

    Notes:
        - Uses matplotlib for plotting.
        - Fills between error bounds if `se` is True.
    """

    # Create the second axis
    ax1 = plt.gca()  # primary axis
    
    # Hide the y-axis labels and ticks for the primary axis
    ax1.yaxis.set_ticklabels([])
    ax1.yaxis.set_ticks([])
    ax2 = ax1.twinx()
    ax2.set_ylim(-ylim_set, ylim_set)

    ax2.yaxis.tick_left()
    cols = {
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

    try:
        spd = DataForPlot(ref_profiles, header)

        ax2.plot(spd.x_axis, [0] * len(spd.x_axis), color="grey", linewidth=0.5)
        header = header.split()[0]
        title = f"Profile for {header}"
        if start is not None and end is not None:
            title += f" from position {start} to {end}"
        ax2.set_title(title)

        if cov:
            if abund:
                spd.convert_to_coverage(abund=True)
            else:
                spd.convert_to_coverage(abund=False)
        if se:
            spd.convert_to_error_bounds()
            spd.flatten(smoothing_window)
            ax2.fill_between(
                spd.x_axis,
                spd.y_flat[0],
                spd.y_flat[2],
                color=cols[spd.srna_len],
                alpha=0.4,
                label=str(spd.srna_len) + " nt",
            )
            ax2.fill_between(
                spd.x_axis,
                spd.y_flat[1],
                spd.y_flat[3],
                color=cols[spd.srna_len],
                alpha=0.4,
            )
        else:
            spd.flatten(smoothing_window)
            first = True
            for i in range(spd.replicates):
                if first:
                    ax2.plot(
                        spd.x_axis,
                        spd.y_flat[i],
                        spd.y_flat[i + spd.replicates],
                        color=cols[spd.srna_len],
                        alpha=0.8,
                        label=str(spd.srna_len) + " nt",
                    )
                    first = False
                else:
                    ax2.plot(
                        spd.x_axis,
                        spd.y_flat[i],
                        spd.y_flat[i + spd.replicates],
                        color=cols[spd.srna_len],
                        alpha=0.8,
                    )
                ax2.fill_between(
                    spd.x_axis,
                    spd.y_flat[i],
                    spd.y_flat[i + spd.replicates],
                    color=cols[spd.srna_len],
                    alpha=0.05,
                )
        if start is not None:
            x_ticks = ax2.get_xticks()
            ax2.set_xticks(
                x_ticks, [int(x + start - padding) for x in x_ticks]
            )  # Adjusting for 30 bp padding
        ax2.yaxis.set_label_coords(-0.1, 0.5)
        handles, labels = ax2.get_legend_handles_labels()
        return handles, labels
    except:
        pass


def comma_separated_ints(value):
    try:
        return [int(v) for v in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Invalid comma-separated integers: '{}'".format(value)
        )


def comma_separated_strings(value):
    return value.split(",")


def align_plot(
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
    """
    Generates an alignment plot for given sequence lengths and saves the plot as SVG.

    Args:
        align_prefix (str): Prefix for alignment files.
        align_lens (list of int): List of sequence lengths to include in the plot.
        header (str or list of str): Header(s) identifying the sequences.
        start (int, optional): Start position for the plot.
        end (int, optional): End position for the plot.
        smoothing_window (int, optional): Smoothing window size. Default is 1.
        cov (bool, optional): Include coverage data in the plot. Default is True.
        abund (bool, optional): Include abundance data in the plot. Default is True.
        se (bool, optional): Include standard error in the plot. Default is True.
        save (bool, optional): Save the plot as an SVG file. Default is True.
        ylim_set (int, optional): Y-axis limit. Default is 0.
        sec_structure (bool, optional): Plot secondary structure. Default is False.
        ref_file (str, optional): Reference file for secondary structure.
        show_seq (bool, optional): Print the selected sequence to the console. Default is False.

    Notes:
        - Assumes the existence of a `RefProfiles` class for loading data.
        - Assumes the existence of a `single_plot` function for plotting individual sequences.
        - Assumes the existence of a `SectStruct` class for handling secondary structures, if applicable.
    """
    ss = None
    if sec_structure:
        ss = SectStruct(ref_file, header, start, end)
        if show_seq:
            print("Selected sequence:")
            print(ss.sequence)

    # Set up the plot
    set_up_plot(ss)

    all_handles = []
    all_labels = []
    if ylim_set ==0:
        ylim_set = get_y_max(align_prefix, align_lens, header, start, end, smoothing_window, cov, abund, se, ylim_set)



    for len in align_lens:
        try:
            file_path = "{0}_{1}.csv".format(align_prefix, len)
            if os.path.isfile(file_path):
                rp = RefProfiles()
                rp.load_single_ref_profiles(file_path, header=header, start=start, end=end)
                if isinstance(header, list):
                    for h in header:
                        ret = single_plot(
                            rp, h, start, end, smoothing_window, cov, abund, se, ylim_set
                        )
                        if ret is not None:
                            handles, labels = ret
                            all_handles.extend(handles)
                            all_labels.extend(labels)
                else:
                    ret = single_plot(
                        rp, header, start, end, smoothing_window, cov, abund, se, ylim_set
                    )
                    if ret is not None:
                        handles, labels = ret
                        all_handles.extend(handles)
                        all_labels.extend(labels)
            else:
                print(f"File {file_path} not found. Skipping.")
        except:
            pass

    if se:
        plt.legend(all_handles, all_labels, loc="upper right")

    if save:
        if isinstance(header, list):
            save_file = align_prefix + "_" + "_".join(header) + ".svg"
        else:
            save_file = align_prefix + "_" + header + ".svg"

        if start is not None and end is not None:
            save_file = save_file.replace(".svg", f"_{start}_{end}.svg")

        plt.savefig(save_file)

    plt.show()

def get_y_max(align_prefix, align_lens, header, start, end, smoothing_window, cov, abund, se, ylim_set):
    for len in align_lens:
        try:
            file_path = "{0}_{1}.csv".format(align_prefix, len)
            if os.path.isfile(file_path):
                rp = RefProfiles()
                rp.load_single_ref_profiles(file_path, header=header, start=start, end=end)

                if isinstance(header, list):
                    for h in header:
                        ylim_set = _parse_alignments(smoothing_window, cov, abund, se, ylim_set, rp, h)
                else:
                    ylim_set = _parse_alignments(smoothing_window, cov, abund, se, ylim_set, rp, header)
            else:
                print(f"File {file_path} not found. Skipping.")
        except:
            pass
    return ylim_set*1.1

def _parse_alignments(smoothing_window, cov, abund, se, ylim_set, rp, h):
    spd = DataForPlot(rp, h)

    if cov:

        if abund:
            spd.convert_to_coverage(abund=True)
        else:
            spd.convert_to_coverage(abund=False)
    
    if se:
        spd.convert_to_error_bounds()
        spd.flatten(smoothing_window)
        ylim_set = max(ylim_set, max(spd.y_flat[2]), abs(min(spd.y_flat[3])))
    else:
        spd.flatten(smoothing_window)
        for i in range(spd.replicates):
            ylim_set = max(ylim_set, max(spd.y_flat[i]), abs(min(spd.y_flat[i + spd.replicates])))
    return ylim_set


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
    parser.add_argument("header", help="Header of plot", type=comma_separated_strings)

    # Optional arguments
    parser.add_argument(
        "-s", "--smoothing_window", help="Smoothing window", type=int, default=1
    )
    parser.add_argument("-c", "--coverage", help="Plot coverage", action="store_true")
    parser.add_argument("-a", "--abundance", help="Plot abundance", action="store_true")
    parser.add_argument("-e", "--error", help="Plot error", action="store_true")
    parser.add_argument(
        "-y", "--ylim", help="Set y-axis limit", type=int, nargs=2, default=0
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

    # Call the align_plot function with the parsed arguments
    align_plot(
        args.align_prefix,
        args.align_lens,
        args.header,
        args.start,
        args.end,
        args.smoothing_window,
        args.coverage,
        args.abundance,
        args.error,
        not args.no_save,
        args.ylim,
        args.sec_structure,
        args.ref_file,
    )


# Note: Add the following lines to trigger the main() function when the script is executed:
if __name__ == "__main__":
    main()
