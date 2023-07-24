import math
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scram2Plot.profilePlot import RefProfiles


def get_bin_data(srp, bin_size):
    """
    This function gets the bin data by taking in two arguments: a single reference profile object srp and a bin_size.
    It generates a bin number based on the length of the reference divided by the bin size.
    It then assigns each alignment to its respective bin based on its position and strand.
    The function returns the populated bins.
    """
    bin_no = math.ceil(srp.ref_len / bin_size)

    # bins is a lists of lists of len bin no., with a f and r list in each list
    bins = np.array(
        [[[0.0] * srp.replicates, [0.0] * srp.replicates] for _ in range(bin_no)]
    )
    for i in srp.all_alignments:
        bin_count = int(i.position / bin_size)
        if i.strand == "+":
            bins[bin_count][0] += i.indv_alignments
        else:
            bins[bin_count][1] += i.indv_alignments
    return bins


def calculate_means_and_ses(bins):
    """
    This function takes in a list of bins and calculates the mean and standard errors
    across all replicates within each bin. The function returns a tuple containing the
    means and standard errors.
    """
    means = np.mean(bins, axis=2)
    ses = np.std(bins, axis=2) / np.sqrt(bins.shape[2])
    return means, ses


def create_dataframe(bins, means, ses):
    """
    This function takes in bins, means, and standard errors, and creates a pandas DataFrame with
    columns for the bin number, list (either 'first' or 'second'), mean, and standard error.
    The DataFrame is then returned.
    """
    data = []
    for i in range(bins.shape[0]):
        data.append([i, "first", means[i, 0], ses[i, 0]])
        data.append([i, "second", -means[i, 1], ses[i, 1]])
    df = pd.DataFrame(data, columns=["bin", "list", "mean", "se"])
    return df


def plot_standard_errors(df, color="gray", alpha=0.3):
    """
    This function plots the standard error region for two lists in a DataFrame (df) with a given color and alpha.
    It returns the DataFrame
    """
    for list_name in ["first", "second"]:
        plt.fill_between(
            df[df["list"] == list_name]["bin"],
            (
                df[df["list"] == list_name]["mean"] - df[df["list"] == list_name]["se"]
            ).values,
            (
                df[df["list"] == list_name]["mean"] + df[df["list"] == list_name]["se"]
            ).values,
            color=color,
            alpha=alpha,
        )
    return df


def create_plot(bins, color="gray", alpha=0.7):
    """
    This function generates a plot for the given bin data.
    It first calculates the means and standard errors, and then creates a DataFrame.
    It then plots the standard error region for the DataFrame.
    The DataFrame is returned.
    """
    means, ses = calculate_means_and_ses(bins)
    df = create_dataframe(bins, means, ses)
    df = plot_standard_errors(df, color, alpha)
    return df


def show_plot(
    dataframes, scale_factor=1, filename=None, y_min=None, y_max=None, cols={}
):
    """
    This function displays the generated plots for a list of dataframes.
    It takes in several parameters including a list of dataframes, a scale factor for the x-axis, a filename (optional),
    minimum and maximum y-axis values (optional), and a dictionary of colors for the legend (optional).
    """
    plt.axhline(0, color="black")
    if y_min is None or y_max is None:  # If y-axis limits are not provided
        y_abs_max = (
            max(max(df["mean"] + df["se"]) for df in dataframes) * 1.1
        )  # Increased limits by 10%
        if y_min is None:
            y_min = -y_abs_max
        if y_max is None:
            y_max = y_abs_max
    plt.ylim(y_min, y_max)
    num_bins = max(df["bin"].max() for df in dataframes) + 1
    ticks = np.linspace(0, num_bins - 1, 10).astype(int)
    scaled_ticks = ticks * scale_factor
    formatted_ticks = [
        format(int(tick), ",") for tick in scaled_ticks
    ]  # Format with commas
    plt.xticks(ticks, formatted_ticks, rotation=45)  # Rotate labels 45 degrees
    plt.xlabel("Reference position (bp)")  # Label for x-axis
    plt.ylabel("Reads per million reads")  # Label for y-axis
    if filename is not None:
        plt.title(filename.split("/")[-1][:-4])
    create_custom_legend(cols)  # create custom legend
    if filename:
        plt.savefig(
            filename, bbox_inches="tight"
        )  # Adjust bounding box to fit rotated labels
    plt.show()


def create_custom_legend(colors_dict, alpha=0.7):
    """
    This function generates a custom legend using a dictionary of colors and an alpha value for transparency.
    The keys in the dictionary are sorted and a patch is created for each key-value pair.
    """
    # Order by converting keys to integers and sorting
    colors_ordered = {k: colors_dict[k] for k in sorted(colors_dict, key=int)}
    patches = [
        mpatches.Patch(color=color, label=str(label), alpha=alpha)
        for label, color in colors_ordered.items()
    ]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc="upper left")


def plot_bin_plot(srps, lens, bin_size, filename=None, y_min=None, y_max=None):
    """
    This function generates a plot of binned data. It takes in several parameters including a list of single
    reference profile objects, a list of lengths, bin size, and optional parameters for filename and y-axis limits.
    """
    plt.figure(figsize=(10, 6))
    cols = {
        24: "darkgreen",
        21: "#984ea3",
        22: "#3374FF",
        18: "#e41a1c",
        19: "#33D4FF",
        20: "#4daf4a",
        23: "#ffff33",
        25: "#a65628",
        26: "#f781bf",
    }
    dfs = []
    pos = 0
    for i in srps:
        dfs.append(create_plot(get_bin_data(i, bin_size), cols[i.srna_len]))
        pos += 1
    show_plot(dfs, bin_size, filename, y_min, y_max, cols)


def load_ref_profiles(file_prefix, lens):
    """
    This function loads reference profiles from a specified file for each length in the lens list.
    It returns a list of reference profiles.
    """
    rp = []
    for i in lens:
        file_name = "{}_{}.csv".format(file_prefix, i)
        rprofile = RefProfiles()
        rprofile.load_single_ref_profiles(file_name)
        rp.append(rprofile)
    return rp


def bin_plot(
    file_prefix, lens, header_prefix, bin_size, out_prefix=None, y_min=None, y_max=None
):
    """
    This function performs the binning and plotting process for a set of reference profiles.
    It loads reference profiles, extracts the single reference profile for each header,
    and then generates a bin plot for each
    """
    rps = load_ref_profiles(file_prefix, lens)
    all_headers = rps[0].single_ref_profiles.keys()
    headers = [i for i in all_headers if i.startswith(header_prefix)]
    for header in headers:
        print("\n" + header)
        srps = []
        for i in rps:
            try:
                srp = i.single_ref_profiles[header]
                srps.append(srp)
            except:
                pass
        if out_prefix is not None:
            filename = "{}_{}.png".format(out_prefix, header.split()[0])
            plot_bin_plot(srps, lens, bin_size, filename, y_min, y_max)
        else:
            plot_bin_plot(srps, lens, bin_size, None, y_min, y_max)


def comma_separated_values(value):
    return [int(i) if i.isdigit() else i for i in value.split(",")]


def main():
    parser = argparse.ArgumentParser(
        description="Plot binned data for long references (e.g. chromosomes)"
    )
    parser.add_argument(
        "align_prefix", type=str, help="Prefix of alignment files"
    )
    parser.add_argument(
        "align_lens",
        type=comma_separated_values,
        help="Comma-separated list of siRNA lengths to plot",
    )
    parser.add_argument(
        "header_prefix", type=str, help="Prefix of reference header/s to plot"
    )
    parser.add_argument("--bin_size", type=int, help="Bin size", default=50000)
    parser.add_argument(
        "--out_prefix", type=str, default=None, help="Plot output file prefix (optional)"
    )
    parser.add_argument(
        "--y_min", type=int, default=None, help="Minimum y-axis value (optional)"
    )
    parser.add_argument(
        "--y_max", type=int, default=None, help="Maximum y-axis value (optional)"
    )

    args = parser.parse_args(sys.argv[1:])

    bin_plot(
        args.align_prefix,
        args.align_lens,
        args.header_prefix,
        args.bin_size,
        args.out_prefix,
        args.y_min,
        args.y_max,
    )


if __name__ == "__main__":
    main()
