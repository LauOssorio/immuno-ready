import pandas as pd
import numpy as np
from src.data_processing.data_loader import load_dataset3_pca

def generate_matrix_for_peptide(peptide, pca_table):
    """
    Converts a peptide sequence into a matrix of PCA-reduced amino acid features.

    Each amino acid in the peptide is mapped to its corresponding PCA vector
    from the provided PCA table.

    Parameters:
    -----------
    peptide : str
        A peptide sequence consisting of amino acid characters.

    pca_table : pd.DataFrame
        A lookup table containing PCA-transformed vectors for each amino acid.

    Returns:
    --------
    np.ndarray
        A 2D NumPy array where each row corresponds to the PCA vector of an amino acid.
    """
    pca_table = load_dataset3_pca()

    return np.stack([pca_table[aa] for aa in peptide])



def generate_matrices_for_dataset(dataset):
    """
    Generates PCA-based feature matrices for a list of peptides in a dataset.

    For each peptide in the dataset, this function applies `generate_matrix_for_peptide`
    to build a matrix representation based on PCA-transformed amino acid properties.

    Parameters:
    -----------
    dataset : pd.DataFrame
        The input DataFrame containing a column 'Epitope - Name' with peptide sequences.

    Returns:
    --------
    list of np.ndarray
        A list where each element is a 2D matrix representing one peptide.
    """
    pca_table = load_dataset3_pca()
    peptides = dataset['Epitope - Name'].tolist()
    return [generate_matrix_for_peptide(p, pca_table) for p in peptides]
