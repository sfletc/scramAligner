import math
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    from .profilePlot import RefProfiles, color_for_length
except ImportError:  # allow running as a plain script: python binPlot.py
    from profilePlot import RefProfiles, color_for_length


def get_bin_data(srp, bin_size):
    """
    This function gets the bin data by taking in two arguments: a single reference profile object srp and a bin_size.
    It generates a bin number based on the length of the reference divided by the bin size.
    It then assigns each alignment to its respective bin based on its position and strand.
    The function returns the populated bins as a (bin_no, 2, replicates) array
    (axis 1 is forward/reverse strand).
    """
    bin_no = math.ceil(srp.ref_len / bin_size)
    bins = np.zeros((bin_no, 2, srp.replicates))
    if srp.positions is None or len(srp.positions) == 0:
        return bins
    # positions are 1-based, so the bin index uses position - 1
    bin_idx = (srp.positions - 1) // bin_size
    strand_idx = np.where(srp.strands, 0, 1)
    np.add.at(bins, (bin_idx, strand_idx), srp.counts)
    return bins


def calculate_means_and_ses(bins):
    """
    This function takes in the binned counts and calculates the mean and standard error
    across all replicates within each bin. The standard error uses the sample standard
    deviation (ddof=1); with a single replicate the standard error is zero.
    The function returns a tuple containing the means and standard errors.
    """
    means = np.mean(bins, axis=2)
    replicates = bins.shape[2]
    if replicates > 1:
        ses = np.std(bins, axis=2, ddof=1) / np.sqrt(replicates)
    else:
        ses = np.zeros_like(means)
    return means, ses


def create_dataframe(bins, means, ses):
    """
    This function takes in bins, means, and standard errors, and creates a pandas DataFrame with
    columns for the bin number, list (either 'first' for the forward strand or 'second' for the
    reverse strand, whose means are negated), mean, and standard error. The DataFrame is then returned.
    """
    bin_no = bins.shape[0]
    df = pd.DataFrame(
        {
            "bin": np.tile(np.arange(bin_no), 2),
            "list": np.repeat(["first", "second"], bin_no),
            "mean": np.concatenate((means[:, 0], -means[:, 1])),
            "se": np.concatenate((ses[:, 0], ses[:, 1])),
        }
    )
    return df


def plot_standard_errors(df, color="gray", alpha=0.3):
    """
    This function plots the mean line and standard error region for both strands in a
    DataFrame (df) with a given color and alpha.
    It returns the DataFrame
    """
    for list_name in ["first", "second"]:
        sub = df[df["list"] == list_name]
        plt.plot(sub["bin"], sub["mean"], color=color, alpha=min(1.0, alpha * 2), linewidth=0.8)
        plt.fill_between(
            sub["bin"],
            (sub["mean"] - sub["se"]).values,
            (sub["mean"] + sub["se"]).values,
            color=color,
            alpha=alpha,
        )
    return df


def create_plot(bins, color="gray", alpha=0.7):
    """
    This function generates a plot for the given bin data.
    It first calculates the means and standard errors, and then creates a DataFrame.
    It then plots the mean and standard error region for the DataFrame.
    The DataFrame is returned.
    """
    means, ses = calculate_means_and_ses(bins)
    df = create_dataframe(bins, means, ses)
    df = plot_standard_errors(df, color, alpha)
    return df


def show_plot(
    dataframes, scale_factor=1, header=None, filename=None, y_min=None, y_max=None, cols=None
):
    """
    This function finalizes and displays the generated plots for a list of dataframes.
    It takes in several parameters including a list of dataframes, a scale factor for the x-axis, a filename (optional),
    minimum and maximum y-axis values (optional), and a dictionary of colors for the legend (optional).
    """

    plt.axhline(0, color="black")

    if y_min is None or y_max is None:  # If y-axis limits are not provided
        # reverse-strand means are negative, so take the largest magnitude of
        # both error bounds on both strands
        y_abs_max = (
            max(
                np.abs(
                    np.concatenate(((df["mean"] + df["se"]), (df["mean"] - df["se"])))
                ).max()
                for df in dataframes
            )
            * 1.1
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
    else:
        plt.title(header)
    create_custom_legend(cols or {})  # create custom legend
    if filename:
        # save before show: interactive backends destroy the figure on show
        plt.savefig(filename, bbox_inches="tight")  # Adjust bounding box to fit rotated labels
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


def plot_bin_plot(srps, bin_size, header, filename=None, y_min=None, y_max=None):
    """
    This function generates a plot of binned data. It takes in several parameters including a list of single
    reference profile objects, a list of lengths, bin size, and optional parameters for filename and y-axis limits.
    """
    fig = plt.figure(figsize=(10, 6))
    try:
        # legend only the lengths actually plotted
        cols = {srp.srna_len: color_for_length(srp.srna_len) for srp in srps}
        dfs = []
        for srp in srps:
            dfs.append(create_plot(get_bin_data(srp, bin_size), cols[srp.srna_len]))
        show_plot(dfs, bin_size, header, filename, y_min, y_max, cols)
    finally:
        plt.close(fig)


def load_ref_profiles(file_prefix, lens):
    """
    This function loads reference profiles from a specified file for each length in the lens list.
    It returns a list of reference profiles.
    """
    rp = []
    for i in lens:
        file_name = "{}_{}.csv".format(file_prefix, i)
        rprofile = RefProfiles()
        rprofile.load_single_ref_profiles(file_name, allow_empty=True)
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
    # Preserve first-seen order while including references unique to later lengths.
    all_headers = dict.fromkeys(header for rp in rps for header in rp.single_ref_profiles)
    headers = [i for i in all_headers if i.startswith(header_prefix)]
    if not headers:
        raise ValueError(
            f"No reference headers starting with '{header_prefix}' found in "
            f"{file_prefix} at lengths {lens}"
        )
    for header in headers:
        srps = [
            i.single_ref_profiles[header] for i in rps if header in i.single_ref_profiles
        ]
        if out_prefix is not None:
            filename = "{}_{}.png".format(out_prefix, header.split()[0])
            plot_bin_plot(srps, bin_size, header, filename, y_min, y_max)
        else:
            plot_bin_plot(srps, bin_size, header, None, y_min, y_max)


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

    try:
        bin_plot(
            args.align_prefix,
            args.align_lens,
            args.header_prefix,
            args.bin_size,
            args.out_prefix,
            args.y_min,
            args.y_max,
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
