import numpy as np
import pandas as pd
from src.config import RAW_DATA_PATH


def load_hla_ligand_atlas(normal_file_name = "hla_2020.12_HLA_aggregated.tsv",
                          metadata_file = "hla_2020.12_HLA_sample_hits.tsv"):
    hla_ligand_atlas_df = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep='\t')
    hla_ligand_atlas_metadata = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep='\t')
    return hla_ligand_atlas_df, hla_ligand_atlas_metadata

def n_indv_per_peptide():
    #TODO: calculate number of donnors in which the peptide was found.
    # #Call the column 'Assay - Number of Subjects Tested'
    pass


def add_relevant_columns(hla_ligand_atlas_df):
    """
    Standardizes the HLA Ligand Atlas dataset to match the format of the cleaned IEDB positive dataset.

    This function:
    - Creates a new DataFrame with required columns matching the cleaned positive dataset
    - Adds constant values for control/healthy sample metadata
    - Maps and filters the MHC restriction class to include only MHC class I and II
    - Returns a cleaned DataFrame ready to be merged with the positive dataset

    Parameters:
    -----------
    hla_ligand_atlas_df : pd.DataFrame
        The raw HLA Ligand Atlas DataFrame, expected to have at least
        'peptide_sequence' and 'hla_class' columns.

    Returns:
    --------
    pd.DataFrame
        A new DataFrame with standardized columns and filtered for MHC class I and II only.
    """
    new_data_frame = pd.DataFrame()

    #Rename columns (peptide-name and MHC restriction class)
    hla_ligand_atlas_df.rename(columns = {'peptide_sequence':'Epitope - Name','hla_class':'MHC Restriction - Class'})

    #Create healthy/control columns
    new_data_frame['Epitope - Name'] = hla_ligand_atlas_df['peptide_sequence']
    new_data_frame['Epitope - Source Organism'] = 'Homo sapiens'
    new_data_frame['Epitope - Species'] = 'Homo sapiens'
    new_data_frame['1st in vivo Process - Process Type'] = 'None'
    new_data_frame['1st in vivo Process - Disease'] = 'Healthy'
    new_data_frame['1st in vivo Process - Disease Stage'] = 'Healthy'
    new_data_frame['Assay - Method'] = 'None'
    new_data_frame['Assay - Response measured'] = 'None'
    new_data_frame['Assay - Qualitative Measure'] = 'Negative'
    new_data_frame['Assay - Response Frequency (%)'] = np.nan

    #Rename values in the hla_class column
    MHC_restriction_map = {
    'HLA-I': 'I',
    'HLA-II': 'II',
    'HLA-I+II': 'non classical'
    }

    new_data_frame['MHC Restriction - Class'] = hla_ligand_atlas_df['hla_class'].map(MHC_restriction_map)

    new_data_frame = new_data_frame[(new_data_frame['MHC Restriction - Class'] == "I")|(new_data_frame['MHC Restriction - Class'] == "II")]


    return new_data_frame


def load_clean_normal():
    df = load_hla_ligand_atlas()
    df = add_relevant_columns(df)
    return df
