import pandas as pd
import numpy as np
from src.data_processing.data_loader import load_dataset3_pca

def generate_matrix_for_peptide(peptide, pca_table):
    pca_table = load_dataset3_pca()

    return np.stack([pca_table[aa] for aa in peptide])



def generate_matrices_for_dataset(dataset):
    pca_table = load_dataset3_pca()
    peptides = dataset['Epitope - Name'].tolist()
    return [generate_matrix_for_peptide(p, pca_table) for p in peptides]
