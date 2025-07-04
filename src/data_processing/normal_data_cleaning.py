import numpy as np
import pandas as pd
from src.data_processing.data_loader import load_raw_hla_ligand_atlas
from src.config import RAW_DATA_PATH



def n_indv_per_peptide():
    peptides_df, metadata_df = load_raw_hla_ligand_atlas()

    full_df = pd.merge(peptides_df, metadata_df,
                       on= ["peptide_sequence_id", "hla_class"])[["peptide_sequence", "donor", "hla_class"]].drop_duplicates()
    full_df["averaged_number_of_subjects"] = full_df.groupby("peptide_sequence")["donor"].transform("nunique")

    full_df = full_df.drop(columns="donor").drop_duplicates()
    return  full_df




def load_clean_normal():
    """
    Standardizes the HLA Ligand Atlas dataset to match the format of the cleaned IEDB positive dataset.

    This function:
    - Creates a new DataFrame with required columns matching the cleaned positive dataset
    - Adds constant values for control/healthy sample metadata
    - Maps and filters the MHC restriction class to include only MHC class I and II
    - Returns a cleaned DataFrame ready to be merged with the positive dataset

    Returns:
    --------
    pd.DataFrame
        A new DataFrame with standardized columns and filtered for MHC class I and II only.
    """
    hla_ligand_atlas_df = n_indv_per_peptide()

    new_data_frame = pd.DataFrame()

    #Rename columns (peptide-name and MHC restriction class)
    hla_ligand_atlas_df.rename(columns = {'peptide_sequence':'Epitope - Name','hla_class':'MHC Restriction - Class'})

    #Create healthy/control columns
    new_data_frame['Epitope - Name'] = hla_ligand_atlas_df['peptide_sequence']
    new_data_frame['averaged_number_of_subjects'] =hla_ligand_atlas_df['averaged_number_of_subjects']
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

    new_data_frame = new_data_frame[(new_data_frame['MHC Restriction - Class'] == "I")]

    new_data_frame.to_csv(RAW_DATA_PATH + "cleaned_hla_ligand_atlas_data.csv", index=False)
    print("Cleaned IEDB data saved to 'cleaned_positive_iedb_data.csv'")

    return new_data_frame


if __name__ == "__main__":
    load_clean_normal()
