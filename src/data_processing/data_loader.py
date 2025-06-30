from src.config import RAW_DATA_PATH
import pandas as pd


def load_hla_ligand_atlas(normal_file_name = "hla_2020.12_HLA_aggregated.tsv",
                          metadata_file = "hla_2020.12_HLA_sample_hits.tsv"):
    hla_ligand_atlas_df = pd.read_csv(RAW_DATA_PATH + normal_file_name, sep='\t')
    hla_ligand_atlas_metadata = pd.read_csv(RAW_DATA_PATH + metadata_file, sep='\t')
    return hla_ligand_atlas_df, hla_ligand_atlas_metadata



if __name__ == "__main__":
    load_hla_ligand_atlas()
