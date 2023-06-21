#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os.path
import csv
import seaborn as sns
import pandas as pd
import numpy as np
import argparse


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

    def load_single_ref_profiles(self, in_file):
        """Loads a scram2 profile file for a single sRNA length

        Args:
            in_file (string): scram2 profile file path
        """
        try:
            with open(in_file, "r") as in_handle:
                reader = csv.reader(in_handle, delimiter=",")
                for row in reader:
                    if row[0] == "Header":
                        continue
                    (
                        header,
                        ref_len,
                        srna,
                        position,
                        strand,
                        times_aligned,
                        *indv_alignments,
                    ) = row

                    # Construct a SingleAlignment object
                    srna = DNA(srna)
                    times_aligned = int(times_aligned)
                    position = int(position)
                    ref_len = int(ref_len)
                    indv_alignments = np.array([float(x) for x in indv_alignments])
                    sa = SingleAlignment(
                        srna, position, strand, times_aligned, indv_alignments
                    )

                    # Add the SingleAlignment to the corresponding SingleRefProfile
                    if header not in self.single_ref_profiles:
                        self.single_ref_profiles[header] = SingleRefProfile()
                        self.single_ref_profiles[header].ref_len = ref_len
                        self.single_ref_profiles[header].replicates = len(
                            indv_alignments
                        )
                        self.single_ref_profiles[header].srna_len = len(srna)
                        self.srna_len = len(srna)
                    self.single_ref_profiles[header].all_alignments.append(sa)

                # Set the number of replicates based on the last read alignment
                if "sa" in locals():
                    self.replicates = len(sa.indv_alignments)
                else:
                    raise ValueError("The input file does not contain any alignments")
        except FileNotFoundError:
            print(f"The input file {in_file} does not exist. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing the input file {in_file}: {e}")


class DataForPlot(object):
    """
    Data for plotting class
    """

    def __init__(self, ref_profiles, header):
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
        """Extracts data from a ref profiles object

        Args:
            ref_profiles (RefProfiles): ref profiles object
            header (string): header of the reference sequence
        """
        for sa in self.ref_profiles.single_ref_profiles[self.header].all_alignments:
            if sa.strand == "+":
                self.fwd[sa.position] = sa.indv_alignments
            else:
                self.rvs[sa.position] = sa.indv_alignments

    def convert_to_coverage(self, abund=True):
        """
        TODO: Fix
        """
        self.fwd = self._coverage_per_strand(self.fwd, abund)
        self.rvs = self._coverage_per_strand(self.rvs, abund)

    def _coverage_per_strand(
        self, old_array, abund
    ):  # TODO: Fix for better use of numpy
        """

        Args:
            old_array (_type_): _description_
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

    def convert_to_error_bounds(self):  # TODO: document and test
        """_summary_"""
        self.fwd = self._error_bounds(self.fwd)
        self.rvs = self._error_bounds(self.rvs)

    def _error_bounds(self, old_array):  # TODO: document and test
        """_summary_"""
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

    def __str__(self):  # TODO: update
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.header,
            self.ref_len,
            self.srna_len,
            self.fwd,
            self.rvs,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):  # TODO: update
        return (
            self.header == other.header
            and self.ref_len == other.ref_len
            and self.srna_len == other.srna_len
            and np.array_equal(self.fwd, other.fwd)
            and np.array_equal(self.rvs, other.rvs)
        )


def align_plot(align_prefix, align_lens, header, smoothing_window=1, cov=True, abund=True, se=True, save=True, ylim_set=(0, 0)):
    file_paths = []
    set_up_plot(ylim_set)
    for i in align_lens:
        file_paths.append("{0}_{1}.csv".format(align_prefix, i))
    rp = RefProfiles()
    for i in file_paths:
        if not os.path.isfile(i):
            pass
        else:
            rp.load_single_ref_profiles(i)
    if isinstance(header, list):
        for h in header:
            single_plot(rp, h, smoothing_window, cov, abund, se)
    else:
        single_plot(rp, header, smoothing_window, cov, abund, se)
    if se:
        plt.legend()
    if save:
        if isinstance(header, list):
            save_file = align_prefix + "_" + "_".join(header) + ".png"
        else:
            save_file = align_prefix + "_" + header + ".png"
        plt.savefig(save_file)
    plt.show()


def set_up_plot(ylim_set):
    """ """
    plt.figure(figsize=(12, 6), dpi=300)
    plt.xlabel("Position")
    plt.ylabel("Abundance")
    plt.title("Abundance Profile")
    if ylim_set != (0, 0):
        plt.ylim(ylim_set[0], ylim_set[1])


def single_plot(ref_profiles, header, smoothing_window, cov, abund, se):
    """ """
    cols = {
        24: "darkgreen",
        21: "red",
        22: "blue",
        18: "#f781bf",
        19: "#a65628",
        20: "#984ea3",
        23: "#999999",
        25: "brown",
        26: "#dede00",
        27: "orange",
        28: "yellow",
    }  # TODO: complete

    try:
        spd = DataForPlot(ref_profiles, header)
        plt.plot(spd.x_axis, [0] * len(spd.x_axis), color="grey", linewidth=0.5)
        if cov:
            if abund:
                spd.convert_to_coverage(abund=True)
            else:
                spd.convert_to_coverage(abund=False)
        if se:
            spd.convert_to_error_bounds()
            spd.flatten(smoothing_window)
            plt.fill_between(
                spd.x_axis,
                spd.y_flat[0],
                spd.y_flat[2],
                color=cols[spd.srna_len],
                alpha=0.4,
                label=str(spd.srna_len) + " nt",
            )
            plt.fill_between(
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
                    plt.plot(
                        spd.x_axis,
                        spd.y_flat[i],
                        spd.y_flat[i + spd.replicates],
                        color=cols[spd.srna_len],
                        alpha=0.8,
                        label=str(spd.srna_len) + " nt",
                    )
                    first = False
                else:
                    plt.plot(
                        spd.x_axis,
                        spd.y_flat[i],
                        spd.y_flat[i + spd.replicates],
                        color=cols[spd.srna_len],
                        alpha=0.8,
                    )
                plt.fill_between(
                    spd.x_axis,
                    spd.y_flat[i],
                    spd.y_flat[i + spd.replicates],
                    color=cols[spd.srna_len],
                    alpha=0.05,
                )
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
    return value.split(',')


# command line interface for profile_plot
def main():
    parser = argparse.ArgumentParser(description="Plot abundance profiles")
    parser.add_argument("align_prefix", help="Prefix of alignment files")
    parser.add_argument(
        "align_lens",
        help="Comma-separated list of siRNA lengths to plot",
        type=comma_separated_ints,
    )
    parser.add_argument("header", help="Header of plot", type=comma_separated_strings)
    parser.add_argument(
        "-s", "--smoothing_window", help="Smoothing window", type=int, default=1
    )
    parser.add_argument("-c", "--coverage", help="Plot coverage", action="store_true")
    parser.add_argument("-a", "--abundance", help="Plot abundance", action="store_true")
    parser.add_argument("-e", "--error", help="Plot error", action="store_true")
    parser.add_argument(
        "-y", "--ylim", help="Set y-axis limit", type=int, nargs=2, default=(0, 0)
    )
    parser.add_argument("-n", "--no_save", help="Do not save plot", action="store_true")
    args = parser.parse_args()
    if args.abundance:
        args.coverage = True
    align_plot(
        args.align_prefix,
        args.align_lens,
        args.header,
        args.smoothing_window,
        args.coverage,
        args.abundance,
        args.error,
        not args.no_save,
        args.ylim,
    )
