import numpy as np
import pandas as pd
from src.data_processing.data_loader import load_raw_hla_ligand_atlas
from src.data_processing.feature_engineering import fill_group_II_status
from src.config import RAW_DATA_PATH



def n_indv_per_peptide():
    """
    Loads raw HLA Ligand Atlas data and computes the number of unique donors
    per peptide sequence, grouped by MHC class.

    Returns:
    --------
    pd.DataFrame
        DataFrame with unique peptide sequences, their MHC class, and
        the number of individuals in which they were observed.
    """
    peptides_df, metadata_df = load_raw_hla_ligand_atlas()

    full_df = pd.merge(peptides_df, metadata_df,
                       on= ["peptide_sequence_id", "hla_class"])[["peptide_sequence", "donor", "hla_class"]].drop_duplicates()
    full_df["averaged_number_positive_subjects_tested"] = full_df.groupby("peptide_sequence")["donor"].transform("nunique")

    full_df = full_df.drop(columns="donor").drop_duplicates()
    return  full_df




def peptide_length(peptide):
    """
    Computes the length of a peptide string.

    Parameters:
    -----------
    peptide : str
        Amino acid sequence of the peptide.

    Returns:
    --------
    int
        Length of the peptide.
    """
    return len(str(peptide))

def drop_large_and_short_sequences (data_frame , min_length, max_length):
    """
    Removes peptide sequences that are shorter or longer than the specified thresholds.

    Parameters:
    -----------
    data_frame : pd.DataFrame
        DataFrame with a column 'Epitope - Name' containing peptide sequences.
    min_length : int
        Minimum allowed peptide length.
    max_length : int
        Maximum allowed peptide length.

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame with peptides in the specified length range.
    """

    # Add peptide length column temporarily
    data_frame.loc[:,'peptide length'] = data_frame['Epitope - Name'].apply(peptide_length)

    # Drop rows with too long sequences
    data_frame = data_frame[data_frame['peptide length'] <= max_length]

    # Drop rows with too short sequences
    data_frame = data_frame[data_frame['peptide length'] >= min_length]

    # Remove temporary column
    data_frame.drop(columns='peptide length', inplace=True)

    return data_frame



def load_clean_normal():
    """
    Cleans and transforms the HLA Ligand Atlas data into a format compatible with
    cleaned IEDB data. This function performs the following:
    - Computes number of individuals per peptide
    - Renames and maps relevant columns
    - Adds standardized metadata fields
    - Filters MHC classes to I and II
    - Filters peptides observed in fewer than 3 individuals
    - Annotates MHC sharing status (I/II)

    Returns:
    --------
    pd.DataFrame
        Cleaned and standardized negative/control dataset.
    """
    hla_ligand_atlas_df = n_indv_per_peptide()

    new_data_frame = pd.DataFrame()

    #Rename columns (peptide-name and MHC restriction class)
    hla_ligand_atlas_df.rename(columns = {'peptide_sequence':'Epitope - Name','hla_class':'MHC Restriction - Class'})

    #Create healthy/control columns
    new_data_frame['Epitope - Name'] = hla_ligand_atlas_df['peptide_sequence']
    new_data_frame['averaged_number_positive_subjects_tested'] = hla_ligand_atlas_df['averaged_number_positive_subjects_tested']
    new_data_frame['1st in vivo Process - Process Type'] = 'None'
    new_data_frame['Assay - Qualitative Measurement'] = 'Negative'


    #Rename values in the hla_class column
    MHC_restriction_map = {
    'HLA-I': 'I',
    'HLA-II': 'II',
    'HLA-I+II': 'non classical'
    }

    new_data_frame['MHC Restriction - Class'] = hla_ligand_atlas_df['hla_class'].map(MHC_restriction_map)



    new_data_frame = new_data_frame[(new_data_frame['MHC Restriction - Class'] == "I")| (new_data_frame['MHC Restriction - Class'] == "II")]

    # Drop peptides found in 1 or 2 individuals
    new_data_frame = new_data_frame[new_data_frame['averaged_number_positive_subjects_tested'] > 2]

    new_data_frame = fill_group_II_status(new_data_frame).drop_duplicates()

    return new_data_frame


if __name__ == "__main__":
    load_clean_normal()
