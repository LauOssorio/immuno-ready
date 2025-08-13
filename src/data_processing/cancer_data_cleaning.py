import numpy as np
import pandas as pd
from src.config import RAW_DATA_PATH
from src.data_processing.data_loader import load_raw_cancer
from src.data_processing.feature_engineering import fill_group_II_status


def select_columns_and_clean_cancer(data_frame):
    """
    Selects relevant columns from the IEDB cancer bench marking dataset.

    """

    list_columns = ['Epitope - Name',
                '1st in vivo Process - Process Type',
                '1st in vivo Process - Disease',
                'Assay - Qualitative Measurement',
                'Assay - Number of Subjects Tested',
                'Assay - Response Frequency (%)',
                'MHC Restriction - Class']

    data_frame = data_frame[list_columns].drop_duplicates()
    # data_frame = data_frame[data_frame["1st in vivo Process - Process Type"] != "No immunization"]
    # data_frame = data_frame[data_frame["1st in vivo Process - Disease"] != "healthy"]

    return data_frame



def fix_weird_peptides (data_frame):
    """
    Cleans peptide sequences by:
    - Removing parts after '+' if present.
    - Removing sequences containing invalid amino acids (U or X).

    Parameters:
        data_frame (pd.DataFrame): DataFrame containing a 'Epitope - Name' column.

    Returns:
        pd.DataFrame: Cleaned DataFrame without malformed peptide entries.
    """

    # Remove the right-hand part for peptides that contain a " + "
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    data_frame = data_frame[data_frame['Epitope - Name'].apply(lambda x: isinstance(x, str) and set(x).issubset(valid_aa))]

    return data_frame



def peptide_length(peptide):
    """
    Computes the length of a given peptide.

    Parameters:
        peptide (str): Amino acid sequence.

    Returns:
        int: Length of the peptide.
    """
    return len(str(peptide))


def drop_large_and_short_sequences (data_frame , min_length, max_length):
    """
    Removes peptides with sequence lengths outside the specified range.

    Parameters:
        data_frame (pd.DataFrame): DataFrame with peptide sequences.
        min_length (int): Minimum allowed peptide length.
        max_length (int): Maximum allowed peptide length.

    Returns:
        pd.DataFrame: Filtered DataFrame with sequences within valid length range.
    """
# Function that drops all peptides with sequence length > 25 (or any chosen max_length) AA

    # Add peptide length column temporarily
    data_frame.loc[:,'peptide length'] = data_frame['Epitope - Name'].apply(peptide_length)

    # Drop rows with too long sequences
    data_frame = data_frame[data_frame['peptide length'] <= max_length]

    # Drop rows with too short sequences
    data_frame = data_frame[data_frame['peptide length'] >= min_length]

    # Remove temporary column
    data_frame.drop(columns='peptide length', inplace=True)

    return data_frame

def average_number_of_individuals(data_frame):
    """
    Computes a weighted average of positive individuals tested per peptide,
    accounting for response frequency.

    Parameters:
        data_frame (pd.DataFrame): DataFrame containing assay response data.

    Returns:
        pd.DataFrame: DataFrame with an added 'averaged_number_positive_subjects_tested' column.
    """
    data_frame["positive_subjects_tested"] = data_frame["Assay - Response Frequency (%)"].fillna(100) * 0.01 * data_frame["Assay - Number of Subjects Tested"]

    data_frame["averaged_number_positive_subjects_tested"] = round((
        data_frame.groupby("Epitope - Name")["positive_subjects_tested"]
        .transform(lambda x: x.sum(min_count=1))).fillna(1))

    # drop duplicated lines
    data_frame = data_frame.drop_duplicates()

    return data_frame


def load_clean_cancer (min_length =8, max_length = 25):
    """
    Loads and processes raw IEDB cancer data by applying the full cleaning pipeline:
    - Removes irrelevant or malformed entries.
    - Filters peptides by length.
    - Calculates average assay statistics.
    - Labels peptides shared between MHC class I and II.

    Parameters:
        min_length (int): Minimum peptide length to retain.
        max_length (int): Maximum peptide length to retain.

    Returns:
        pd.DataFrame: Cleaned and annotated IEDB dataset.
    """
    # Loading the IEDB data
    data_frame = load_raw_cancer()

    # Remove unnecessary columns and filter the data
    data_frame = select_columns_and_clean_cancer(data_frame)

    # Remove / fix weird peptides
    data_frame = fix_weird_peptides(data_frame)

    # Remove unusually long sequences (default is >25 but any max_length will work)
    data_frame = drop_large_and_short_sequences(data_frame , min_length, max_length)

    # calculate the averaged number of individuals used in the assays per peptide
    data_frame = average_number_of_individuals(data_frame)

    data_frame = data_frame.dropna(subset=['MHC Restriction - Class'])

    # add mhc group status for peptides that are found in MHC I and II
    data_frame = fill_group_II_status(data_frame)

    data_frame = data_frame[data_frame["1st in vivo Process - Process Type"] == 'Occurrence of cancer']

    data_frame = data_frame.drop(columns=[
        'Assay - Number of Subjects Tested',
                'Assay - Response Frequency (%)',
                '1st in vivo Process - Disease',
                'positive_subjects_tested'
    ]).drop_duplicates()
    data_frame = data_frame[(data_frame['MHC Restriction - Class'] == "I")| (data_frame['MHC Restriction - Class'] == "II")]



    return data_frame
