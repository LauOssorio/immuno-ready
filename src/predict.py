import joblib
import numpy as np
import re

from tensorflow.keras.models import load_model
import pandas as pd
from src.config import SAVED_MODELS_PATH

from src.data_processing.sequence_tokenizer import AA_index_tokenizer


def predict_new_samples_cnn_multimodal_classificator(X_new, X_cat_new, threshold=0.4):
    """
    Predicts immunogenic classification outcomes for new peptide samples using a pretrained
    multimodal CNN model with fixed architecture.

    This function:
    1. Filters out peptide sequences that are not within the valid amino acid range (9–20 residues, using standard AA codes).
    2. Tokenizes filtered peptide sequences using the AA_index-based tokenizer.
    3. Preprocesses the associated categorical metadata using one-hot encoding.
    4. Loads a trained Keras model from file.
    5. Predicts classification probabilities and applies a threshold to generate binary labels.

    Args:
        model_path (str): Path to the trained Keras model cnn_multimodal_class.h5.
        X_new (pd.DataFrame): DataFrame containing peptide sequences in the 'Epitope Name' column.
        X_cat_new (pd.DataFrame): DataFrame containing associated categorical metadata
                                  (must include 'MHC Restriction - Class', 'mhc_status', and be aligned with X_new).
        threshold (float): Probability threshold for binary classification (default is 0.4).

    Returns:
        pd.DataFrame: DataFrame with two columns:
            - 'target_prob': Predicted probabilities of immunogenicity.
            - 'target_strength': Binary classification labels (1 = positive, 0 = negative).

    Raises:
        ValueError: If required columns are missing from input DataFrames.

    """
    model_path = SAVED_MODELS_PATH + "cnn_multimodal_class.h5"
    valid_aa_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]{9,20}$')

    df_filtered = X_new[X_new['Epitope Name'].apply(lambda x: bool(valid_aa_pattern.match(str(x))))]
    X_new_tokenized = AA_index_tokenizer(df_filtered)
    categorical_features = ['MHC Restriction - Class', 'mhc_status']
    X_cat_new_processed= pd.get_dummies(X_cat_new, columns=categorical_features).drop(columns = ['Epitope - Name', 'peptide_strength','peptide_base'])

    model = load_model(model_path)
    pred_probs = model.predict([X_new_tokenized, X_cat_new_processed])
    y_pred_prob = pred_probs.flatten()
    y_pred_label = (y_pred_prob > threshold).astype(int)

    return pd.DataFrame({
        'target_prob': y_pred_prob,
        'target_strength': y_pred_label
    })
