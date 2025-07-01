import numpy as np
import pandas as pd
from src.config import RAW_DATA_PATH
from src.data_processing.data_loader import load_iedb


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
                'Epitope - Source Organism',
                'Epitope - Species',
                '1st in vivo Process - Process Type',
                '1st in vivo Process - Disease',
                '1st in vivo Process - Disease Stage',
                'Assay - Method',
                'Assay - Response measured',
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

# TODO: the response frequency is not consisted with the sequences becasue the same sequence
# was tested in different assays - do something with this so one peptide will have one
# frequency value. Relate to the averaged number of subjects.


def peptide_length(peptide):
# Function to calculate peptide length
        return len(str(peptide))


def drop_large_and_short_sequences (data_frame , min_length, max_length):
# Function that drops all peptides with sequence length > 25 (or any chosen max_length) AA
# FYI: choosing 25 removes c.7k rows (after removal of weird values)

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
    data_frame["averaged_number_subjects_tested"] = (
        data_frame.groupby("Epitope - Name")["Assay - Number of Subjects Tested"]
        .transform(lambda x: x.sum(min_count=1)))

    data_frame.drop(columns=["Assay - Number of Subjects Tested"], inplace=True)

    # drop duplicateed lines
    data_frame = data_frame.drop_duplicates()

    return data_frame


def load_clean_iedb (min_length =8, max_length = 25):
    # Final cleaning function
    # Loading the IEDB data
    data_frame = load_iedb()

    # Remove unnecessary columns and filter the data
    data_frame = select_columns_and_clean_iedb(data_frame)

    # Remove / fix weird peptides
    data_frame = fix_weird_peptides (data_frame)

    # Remove overlap between neg and pos (only negatives are deleted)
    data_frame = remove_overlap (data_frame)

    # Remove unusually long sequences (default is >25 but any max_length will work)
    data_frame = drop_large_and_short_sequences (data_frame , min_length, max_length)

    # calculate the averaged number of individuals used in the assays per peptide
    data_frame = average_number_of_individuals(data_frame)

    return data_frame


if __name__ == "__main__":
    # Load the IEDB data and clean it
    cleaned_data = load_clean_iedb()

    # Save the cleaned data to a CSV file
    cleaned_data.to_csv(RAW_DATA_PATH + "cleaned_positive_iedb_data.csv", index=False)
    print("Cleaned IEDB data saved to 'cleaned_positive_iedb_data.csv'")
