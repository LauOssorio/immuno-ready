from src.config import RAW_DATA_PATH
from src.config import PROCESSED_DATA_PATH

import pandas as pd


def load_raw_hla_ligand_atlas(normal_file_name = "hla_2020.12_HLA_aggregated.tsv",
                          metadata_file = "hla_2020.12_HLA_sample_hits.tsv"):
    """
    Load HLA Ligand Atlas data and associated metadata.
    """
    hla_ligand_atlas_df = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep='\t')
    hla_ligand_atlas_metadata = pd.read_csv(RAW_DATA_PATH + metadata_file, sep='\t')
    return hla_ligand_atlas_df, hla_ligand_atlas_metadata

def load_raw_iedb(positive_file_name = "tcell_table_export_1751306060.csv"):
    """
    Load IEDB T-cell dataset from CSV.
    """
    iedb_df = pd.read_csv(RAW_DATA_PATH + positive_file_name, sep=',')
    return iedb_df


def load_raw_cancer(cancer_file_name = "benchmark_cancer_positive_negative_tcell_table_export_1753627681.csv"):
    """
    Load IEDB cancer T-cell dataset (positive and negative results) from CSV.
    """
    cancer_df = pd.read_csv(RAW_DATA_PATH + cancer_file_name, sep=',', low_memory=False)
    return cancer_df

def load_dataset3_pca(file_name="dataset3_pca.csv"):
    """
    Load the PCA-reduced AA index dataset developed by Ben Galet, PhD, from a CSV file.
    """
    pca_df = pd.read_csv(PROCESSED_DATA_PATH + file_name, sep=',')
    return pca_df
