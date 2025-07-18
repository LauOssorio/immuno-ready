from src.config import RAW_DATA_PATH
from src.config import PROCESSED_DATA_PATH

import pandas as pd


def load_raw_hla_ligand_atlas(normal_file_name = "hla_2020.12_HLA_aggregated.tsv",
                          metadata_file = "hla_2020.12_HLA_sample_hits.tsv"):
    hla_ligand_atlas_df = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep='\t')
    hla_ligand_atlas_metadata = pd.read_csv(RAW_DATA_PATH + metadata_file, sep='\t')
    return hla_ligand_atlas_df, hla_ligand_atlas_metadata

def load_raw_iedb(normal_file_name = "tcell_table_export_1751306060.csv"):
    iedb_df = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep=',')
    return iedb_df

def load_dataset3_pca(file_name="dataset3_pca.csv"):
    """
    Load the PCA-reduced dataset from a CSV file.
    """
    pca_df = pd.read_csv(PROCESSED_DATA_PATH + file_name, sep=',')
    return pca_df

if __name__ == "__main__":
    def load_all_rawdata():
        df_ligand = load_raw_hla_ligand_atlas()
        df_iedb = load_raw_iedb()
        return df_ligand, df_iedb
    load_all_rawdata()
