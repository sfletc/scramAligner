import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scram2Plot.profilePlot as pp


warnings.filterwarnings("ignore", message=".*defaulting to pandas implementation*")
warnings.filterwarnings("ignore", message=".*This may take some time*")


class ScramCSV:
    """
    Class for handling scram aligner CSV files and generating sequence logos.
    
    Attributes:
        df (DataFrame): The DataFrame containing the loaded data.
        engine (str): The DataFrame engine to use ('pandas' or 'modin').
        pd (module): The DataFrame library in use (either pandas or modin.pandas).
    """
        
    def __init__(self, engine="pandas"):
        """
        Initializes a new instance of the ScramCSV class.
        
        Args:
            engine (str): The DataFrame engine to use ('pandas' or 'modin').
        """
        self.df = None
        self.engine = engine

        if self.engine == "modin":
            import modin.pandas as pd
        else:
            import pandas as pd
        self.pd = pd
        self.srna_len=0

    def load_csv(self, file_path):
        """
        Loads a CSV file into a DataFrame.
        
        Args:
            file_path (str): The path to the CSV file to load.
        """
        self.df = self.pd.read_csv(file_path)
        self.srna_len=int(file_path.split("_")[-1].split(".csv")[0])

    def subset_data(
        self,
        header=None,
        start_position=None,
        end_position=None,
        strand=None,
        times_aligned_min=None,
        times_aligned_max=None,
        first_nucleotide=None,
        min_count=None,
        max_count=None,
    ):
        """
        Filters the data based on various conditions.
        
        Args:
            header (str, optional): The header to filter by.
            start_position (int, optional): The start position to filter by.
            end_position (int, optional): The end position to filter by.
            strand (str, optional): The strand to filter by.
            times_aligned_min (int, optional): The minimum times aligned to filter by.
            times_aligned_max (int, optional): The maximum times aligned to filter by.
            first_nucleotide (str, optional): The first nucleotide to filter by.
            min_count (int, optional): The minimum count to filter by.
            max_count (int, optional): The maximum count to filter by.
        
        Returns:
            ScramCSV: A new ScramCSV instance containing the filtered data.
        """
        if self.df is None:
            raise ValueError("No CSV data loaded. Use load_csv to load a CSV file.")

        # index the mask like the DataFrame so filtering a filtered instance
        # (whose index is no longer a fresh RangeIndex) still aligns correctly
        mask = self.pd.Series(True, index=self.df.index)

        # Apply each filter on the mask
        if header is not None:
            mask &= self.df["Header"] == header
        if start_position is not None:
            mask &= self.df["Position"] >= start_position
        if end_position is not None:
            mask &= self.df["Position"] <= end_position
        if strand is not None:
            mask &= self.df["Strand"] == strand
        if times_aligned_min is not None:
            mask &= self.df["Times aligned"] >= times_aligned_min
        if times_aligned_max is not None:
            mask &= self.df["Times aligned"] <= times_aligned_max
        if first_nucleotide is not None:
            mask &= self.df["sRNA"].str.startswith(first_nucleotide)

        # Filter DataFrame using the final mask
        filtered_df = self.df[mask].copy()

        # Perform other operations such as computing the mean count
        if min_count is not None or max_count is not None:
            count_columns = filtered_df.columns[
                filtered_df.columns.get_loc("Times aligned") + 1 :
            ]
            if len(count_columns) == 0:
                raise ValueError(
                    "No read-count columns found after 'Times aligned'; count "
                    "filtering requires the aligner's per-replicate CSV format"
                )
            filtered_df["Mean_Count"] = filtered_df[count_columns].mean(axis=1)

            if min_count is not None:
                filtered_df = filtered_df[filtered_df["Mean_Count"] >= min_count]
            if max_count is not None:
                filtered_df = filtered_df[filtered_df["Mean_Count"] <= max_count]
            filtered_df = filtered_df.drop(columns=["Mean_Count"])

        # Create a new instance of the class with the filtered DataFrame
        filtered_instance = self.__class__(engine=self.engine)
        filtered_instance.df = filtered_df
        filtered_instance.srna_len = self.srna_len
        return filtered_instance

    def remove_duplicates(self):
        """
        Removes duplicate rows based on the 'sRNA' column.
        """
        if self.df is None:
            raise ValueError(
                "No data to remove duplicates from. Use load_csv or subset_data to load or filter data."
            )

        self.df = self.df.drop_duplicates(subset=["sRNA"])

    def generate_sequence_logo(self, save_path=None, logo_type="information"):
        """
        Generates and displays a sequence logo from the data.
        
        Args:
            save_path (str, optional): The file path to save the generated logo to.
        """
        if self.df is None:
            raise ValueError(
                "No data to generate a sequence logo from. Use load_csv or subset_data to load or filter data."
            )

        import logomaker as lm  # deferred: heavy optional dependency

        # Extract 'sRNA' sequences from the DataFrame
        seqs = self.df["sRNA"].tolist()
        counts_mat = lm.alignment_to_matrix(seqs, to_type="counts")
        if logo_type == "information":
            info_mat = lm.transform_matrix(
                counts_mat, from_type="counts", to_type="information"
            )

            logo = lm.Logo(info_mat)
        else:
            logo = lm.Logo(counts_mat)

        # Ensure that x-axis ticks are integers starting at 1
        xticks = np.arange(self.srna_len)  # Create an array of tick positions corresponding to the length of your data
        xticks_labels = xticks + 1  # Add 1 to each tick position for the label

        plt.xticks(xticks, labels=xticks_labels.astype(int))  # Set tick positions and labels
        # Display the logo
        # Add title using matplotlib
        plt.title("{1} nt sRNAs ({0} input sequences)".format(len(self.df),self.srna_len))

        # Show the plot
        plt.show()

        # Save the logo if save_path is specified
        if save_path is not None:
            logo.savefig(save_path)

    def save_to_csv(self, file_path):
        """
        Saves the data to a CSV file.
        
        Args:
            file_path (str): The file path to save the data to.
        """
        if self.df is None:
            raise ValueError(
                "No data to save. Use load_csv or subset_data to load or filter data."
            )

        self.df.to_csv(file_path, index=False)

    def calculate_min_free_energy(self, ref_file, down=100, up=100):
        """
        Calculate the minimum free energy of the folded window around each sRNA.

        The reference FASTA is indexed once (not re-parsed per row), windows are
        clamped to the sequence bounds, and identical windows are folded only once.

        Args:
            ref_file (str): Path to the reference FASTA file.
            down (int): Bases to include upstream of the alignment position.
            up (int): Bases to include downstream of the alignment position.
        """
        if self.df is None:
            raise ValueError(
                "No data to operate on. Use load_csv or subset_data to load or filter data."
            )

        from Bio import SeqIO
        from ViennaRNA import RNA

        ref_index = SeqIO.index(ref_file, "fasta")
        sequences = {}
        fold_cache = {}
        energies = []
        for header, position in zip(self.df["Header"], self.df["Position"]):
            if header not in sequences:
                if header not in ref_index:
                    raise ValueError(
                        f"Header {header} not found in reference file {ref_file}"
                    )
                sequences[header] = str(ref_index[header].seq)
            seq = sequences[header]
            window = seq[max(0, position - down) : min(len(seq), position + up)]
            if window not in fold_cache:
                fold_cache[window] = RNA.fold(window)[1]
            energies.append(fold_cache[window])
        self.df["Min_Free_Energy"] = energies


