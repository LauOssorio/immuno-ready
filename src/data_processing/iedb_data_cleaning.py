import numpy as np
import pandas as pd
from src.config import RAW_DATA_PATH
from src.data_processing.data_loader import load_raw_iedb
from src.data_processing.feature_engineering import fill_group_II_status


def select_columns_and_clean_iedb(data_frame):
    """
    Filters and cleans an IEDB dataset by selecting relevant columns and removing
    observations that are not associated with immunization or disease.

    This function performs the following steps:
    1. Selects a predefined set of biologically relevant columns related to epitopes,
       assay conditions, and MHC restrictions.
    2. Drops duplicate rows to avoid redundant entries.
    3. Removes rows corresponding to observations from healthy tissues or cases with
       no immunization, as they are not informative for immunogenicity analysis.

    Parameters:
    ----------
    data_frame : pandas.DataFrame
        The input DataFrame containing IEDB data.

    Returns:
    -------
    pandas.DataFrame
        A cleaned DataFrame containing only relevant columns and filtered observations.
    """

    list_columns = ['Epitope - Name',
                '1st in vivo Process - Process Type',
                '1st in vivo Process - Disease',
                'Assay - Qualitative Measurement',
                'Assay - Number of Subjects Tested',
                'Assay - Response Frequency (%)',
                'MHC Restriction - Class']

    data_frame = data_frame[list_columns].drop_duplicates()
    data_frame = data_frame[data_frame["1st in vivo Process - Process Type"] != "No immunization"]
    data_frame = data_frame[data_frame["1st in vivo Process - Disease"] != "healthy"]

    return data_frame



def fix_weird_peptides (data_frame):
# Remove weird peptides values

    # Remove the right-hand part for peptides that contain a " + "
    data_frame.loc[:,'Epitope - Name'] = data_frame['Epitope - Name'].str.split('+').str[0].str.strip()

    # Boolean mask to filter out sequences that contain "U"
    u_mask = data_frame['Epitope - Name'].str.contains('U')
    data_frame = data_frame[~u_mask]

    # Boolean mask to filter out sequences that contain "X"
    x_mask = u_mask = data_frame['Epitope - Name'].str.contains('X')
    data_frame = data_frame[~x_mask]

    return data_frame



def peptide_length(peptide):
# Function to calculate peptide length
        return len(str(peptide))


def drop_large_and_short_sequences (data_frame , min_length, max_length):
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
# Function to calculate the average number of individuals tested per peptide
# correcting the number of individuals by the % of response frequency.
    data_frame["positive_subjects_tested"] = data_frame["Assay - Response Frequency (%)"].fillna(100) * 0.01 * data_frame["Assay - Number of Subjects Tested"]

    data_frame["averaged_number_positive_subjects_tested"] = round((
        data_frame.groupby("Epitope - Name")["positive_subjects_tested"]
        .transform(lambda x: x.sum(min_count=1))).fillna(1))

    # drop duplicated lines
    data_frame = data_frame.drop_duplicates()

    return data_frame


def load_clean_iedb (min_length =8, max_length = 25):
    # Final cleaning function
    # Loading the IEDB data
    data_frame = load_raw_iedb()

    # Remove unnecessary columns and filter the data
    data_frame = select_columns_and_clean_iedb(data_frame)

    # Remove / fix weird peptides
    data_frame = fix_weird_peptides (data_frame)

    # Remove unusually long sequences (default is >25 but any max_length will work)
    data_frame = drop_large_and_short_sequences (data_frame , min_length, max_length)

    # calculate the averaged number of individuals used in the assays per peptide
    data_frame = average_number_of_individuals(data_frame)

    data_frame = data_frame.dropna(subset=['MHC Restriction - Class'])

    # add mhc group status for peptides that are found in MHC I and II
    data_frame = fill_group_II_status(data_frame)

    data_frame = data_frame.drop(columns=[
        'Assay - Number of Subjects Tested',
                'Assay - Response Frequency (%)',
                '1st in vivo Process - Disease',
                'positive_subjects_tested'
    ]).drop_duplicates()
    data_frame = data_frame[(data_frame['MHC Restriction - Class'] == "I")| (data_frame['MHC Restriction - Class'] == "II")]



    return data_frame


if __name__ == "__main__":
    # Load the IEDB data and clean it
    load_clean_iedb()
