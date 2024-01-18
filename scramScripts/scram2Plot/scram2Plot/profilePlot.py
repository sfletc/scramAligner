import numpy as np
import matplotlib.pyplot as plt
import os.path
import csv
import argparse
from Bio import SeqIO
from ViennaRNA import RNA


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
        Sets up a matplotlib plot for displaying biological sequence data.
        """
        # Create a new figure with specified size and DPI
        plt.figure(figsize=(12, 12), dpi=600)

        # Set the x-axis label
        plt.xlabel("Position")

        # Make the x and y dimensions equal; TODO: Check if this is needed
        plt.axis("equal")

        # If a secondary structure object is provided, perform additional setup
        if self.ss is not None:
            self.sec_struct_setup()  # Note: `sec_struct_setup` needs to be defined as another method in this class

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

        self.set_up_plot()

        if self.ylim_set == 0:
            self.ylim_set = self.get_y_max()

        for len in self.align_lens:
            self.process_single_alignment(len)

        if self.se:
            plt.legend(self.all_handles, self.all_labels, loc="upper right")

        if self.save:
            self.save_plot()

        plt.show()

    def get_y_max(self):
        """
        Get the maximum value for the y-axis based on the data.

        Returns:
            float: The maximum y-axis value multiplied by 1.1 for padding.
        """
        for len in self.align_lens:
            self.process_single_alignment(len, get_max=True)
        return self.ylim_set * 1.1

    def process_single_alignment(self, srna_len, get_max=False):
        """
        Processes a single alignment file.

        Args:
            len (int): The length of the alignment.
            get_max (bool, optional): Whether to only get the maximum y-axis value. Defaults to False.
        """
        try:
            file_path = f"{self.align_prefix}_{srna_len}.csv"
            if os.path.isfile(file_path):
                rp = RefProfiles()
                rp.load_single_ref_profiles(
                    file_path, header=self.header, start=self.start, end=self.end
                )

                if isinstance(self.header, list):
                    for h in self.header:
                        self.process_header(rp, h, get_max)
                else:
                    self.process_header(rp, self.header, get_max)
            else:
                print(f"File {file_path} not found. Skipping.")
        except Exception as e:
            print(f"An error occurred: {e}")

    def process_header(self, rp, h, get_max):
        """
        Processes the header information for a given alignment.

        Args:
            rp (RefProfiles): An instance of the RefProfiles class.
            h (str): The header of the sequence.
            get_max (bool): Whether to only get the maximum y-axis value.
        """
        spd = DataForPlot(rp, h)
        self._parse_alignments(spd, h, get_max)

    def _parse_alignments(self, spd, h, get_max):
        """
        Private method for parsing alignments.

        Args:
            spd (DataForPlot): An instance of the DataForPlot class.
            h (str): The header of the sequence.
            get_max (bool): Whether to only get the maximum y-axis value.
        """
        if self.cov:
            spd.convert_to_coverage(abund=self.abund)

        if self.se:
            spd.convert_to_error_bounds()
            spd.flatten(self.smoothing_window)
            if get_max:
                self.ylim_set = max(
                    self.ylim_set, max(spd.y_flat[2]), abs(min(spd.y_flat[3]))
                )
            else:
                self.single_plot(spd, h)  # Call single_plot only if get_max is False
        else:
            spd.flatten(self.smoothing_window)
            if get_max:
                for i in range(spd.replicates):
                    self.ylim_set = max(
                        self.ylim_set,
                        max(spd.y_flat[i]),
                        abs(min(spd.y_flat[i + spd.replicates])),
                    )
            else:
                self.single_plot(spd, h)  # Call single_plot only if get_max is False

    def single_plot(self, spd, header):
        """
        Plots the data for a single sequence based on its reference profiles.

        Args:
            spd (DataForPlot): An instance of the DataForPlot class.
            header (str): The header of the sequence.
        """
        ax1 = plt.gca()  # primary axis
        ax2 = ax1.twinx()  # secondary axis

        self.configure_plot_axes(ax1, ax2)
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
            ax2.plot(spd.x_axis, [0] * len(spd.x_axis), color="grey", linewidth=0.5)
            self.set_plot_title(ax2, header)

            if self.cov:
                spd.convert_to_coverage(abund=self.abund)

            if self.se:
                self.plot_with_error_bounds(ax2, spd, cols)
            else:
                self.plot_reps_individually(ax2, spd, cols)

            if self.start is not None:
                x_ticks = ax2.get_xticks()
                ax2.set_xticks(x_ticks, [int(x + self.start) for x in x_ticks])

            ax2.yaxis.set_label_coords(-0.1, 0.5)
            handles, labels = ax2.get_legend_handles_labels()
            self.all_handles.extend(handles)
            self.all_labels.extend(labels)
            return handles, labels

        except Exception as e:
            print(f"An error occurred while plotting: {e}")

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

    def plot_with_error_bounds(self, ax2, spd, cols):
        """
        Plots the data with error bounds.

        Args:
            ax2: The secondary axis.
            spd (DataForPlot): An instance of the DataForPlot class.
            cols (dict): Dictionary of colors.
        """
        spd.convert_to_error_bounds()
        spd.flatten(self.smoothing_window)
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

    def plot_reps_individually(self, ax2, spd, cols):
        """
        Plots the data without error bounds.

        Args:
            ax2: The secondary axis.
            spd (DataForPlot): An instance of the DataForPlot class.
            cols (dict): Dictionary of colors.
        """
        spd.flatten(self.smoothing_window)
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
    Saves the generated plot as an SVG file.
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
    plotter.generate_plot()


# Note: Add the following lines to trigger the main() function when the script is executed:
if __name__ == "__main__":
    main()
