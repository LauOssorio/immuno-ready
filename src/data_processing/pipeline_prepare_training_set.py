import pandas as pd
import numpy as np
import joblib
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.sequence import pad_sequences


from src.config import PROCESSED_DATA_PATH


from src.data_processing.iedb_data_cleaning import load_clean_iedb
from src.data_processing.normal_data_cleaning import load_clean_normal
from src.data_processing.target_engineering import create_target_features
from src.data_processing.sequence_tokenizer import generate_matrices_for_dataset

def prepare_training_set():
    """
    Prepares the training dataset by merging and processing positive and negative examples,
    engineering features, and saving preprocessed matrices for future modeling.

    This function performs the following steps:
    1. Loads and cleans the IEDB (positive) and HLA Ligand Atlas (negative/control) datasets.
    2. Combines both datasets into a unified DataFrame.
    3. Generates binary target labels using `create_target_features`.
    4. Encodes categorical features (`MHC Restriction - Class`, `mhc_status`) using one-hot encoding.
    5. Generates peptide-level PCA-transformed amino acid features via `generate_matrices_for_dataset`.
    6. Extracts encoded categorical features from the DataFrame.
    7. Separates the target variable (`target_strength`) and sample weights
       (`averaged_number_positive_subjects_tested`).
    8. Saves all processed components as `.joblib` files for reuse during model training.

    Outputs:
    --------
    - X_pca_aa.joblib: numpy array of PCA-transformed amino acid features
    - X_categorical.joblib: numpy array of one-hot encoded categorical features
    - Y.joblib: target labels
    - sample_weights.joblib: sample weights for each observation
    """

    normals  = load_clean_normal()
    iedb_df = load_clean_iedb()
    # Combine the IEDB and normal dataframes
    full_df= pd.concat([normals, iedb_df], ignore_index=True)
    target_df = create_target_features(full_df)

    # Define the categorical features to be one-hot encoded
    categorical_features = ['MHC Restriction - Class', 'mhc_status']

    # One-hot encode the categorical features and append them to the dataframe
    target_encoded = pd.get_dummies(target_df, columns=categorical_features)

    # Generate amino acid PCA embeddings for each peptide sequence (as a list of matrices)
    X_pca_aa = generate_matrices_for_dataset(target_encoded)

    # Extract the last 4 columns which correspond to the one-hot encoded categorical variables
    X_categorical = target_encoded.iloc[:, -4:]  # Assumes exactly 2 categorical features with 2 levels each

    # Extract the target variable representing immunogenicity strength
    Y = target_encoded['target_strength']

    # Extract the raw sample weights (number of individuals supporting the observation)
    sample_weights = target_encoded['averaged_number_positive_subjects_tested']

    # Scale the sample weights to a fixed range (e.g., [0.1, 0.2]) to avoid large disparities in loss contribution
    scaler = MinMaxScaler(feature_range=(0.1, 0.2))
    array_weights = sample_weights.to_numpy().reshape(-1, 1)
    scaled_sample_weights = scaler.fit_transform(array_weights).flatten()

    # Save the processed datasets and scaled weights for future reuse
    joblib.dump(X_pca_aa, PROCESSED_DATA_PATH + 'X_pca_aa.joblib', compress=3)
    joblib.dump(X_categorical, PROCESSED_DATA_PATH + 'X_categorical.joblib', compress=3)
    joblib.dump(Y, PROCESSED_DATA_PATH + 'Y.joblib', compress=3)
    joblib.dump(scaled_sample_weights, PROCESSED_DATA_PATH + 'scaled_sample_weights.joblib', compress=3)





def load_or_create_training_data():
    """
    Loads preprocessed training data from disk if available; otherwise, generates it.

    This function checks for the existence of four joblib files:
    - X_pca_aa (PCA-reduced amino acid features)
    - X_categorical (encoded categorical features)
    - Y (target values)
    - scaled_sample_weights (normalized sample weights)

    If all files exist, they are loaded from disk. If any are missing,
    `prepare_training_set()` is run to generate and save the files before loading.

    Returns:
        tuple: (X_pca_aa, X_categorical, Y, scaled_sample_weights)
    """
    files = {
        "X_pca_aa": PROCESSED_DATA_PATH + "X_pca_aa.joblib",
        "X_categorical": PROCESSED_DATA_PATH + "X_categorical.joblib",
        "Y": PROCESSED_DATA_PATH + "Y.joblib",
        "sample_weights": PROCESSED_DATA_PATH + "scaled_sample_weights.joblib"
    }

    # Check if all files exist
    if all(os.path.exists(path) for path in files.values()):
        print("✅ All joblib files found. Loading...")
        X_pca_aa = joblib.load(files["X_pca_aa"])
        X_categorical = joblib.load(files["X_categorical"])
        Y = joblib.load(files["Y"])
        sample_weights = joblib.load(files["sample_weights"])
    else:
        print("⚠️ One or more files missing. Running prepare_training_set()...")
        prepare_training_set()
        # Now load them
        X_pca_aa = joblib.load(files["X_pca_aa"])
        X_categorical = joblib.load(files["X_categorical"])
        Y = joblib.load(files["Y"])
        sample_weights = joblib.load(files["sample_weights"])

    X_pca_aa_pad = pad_sequences(X_pca_aa, dtype='float32', padding='post')

    return X_pca_aa_pad, X_categorical, Y, sample_weights




def separate_training_test():
    """
    Loads or creates the training data and splits it into training, validation, and test sets (60/20/20).

    Returns:
        tuple: X_pca_train, X_pca_val, X_pca_test,
               X_cat_train, X_cat_val, X_cat_test,
               y_train, y_val, y_test,
               w_train, w_val, w_test
    """
    X_pca_aa_pad, X_categorical, Y, scaled_sample_weights = load_or_create_training_data()

    # First split: Train vs Temp (temp will be split into val and test)
    X_pca_train, X_pca_temp, X_cat_train, X_cat_temp, y_train, y_temp, w_train, w_temp = train_test_split(
        X_pca_aa_pad, X_categorical, Y, scaled_sample_weights,
        test_size=0.4, random_state=42
    )

    # Second split: Validation and Test (from temp)
    X_pca_val, X_pca_test, X_cat_val, X_cat_test, y_val, y_test, w_val, w_test = train_test_split(
        X_pca_temp, X_cat_temp, y_temp, w_temp,
        test_size=0.5, random_state=42
    )

    return (
        X_pca_train, X_pca_val, X_pca_test,
        X_cat_train, X_cat_val, X_cat_test,
        y_train, y_val, y_test,
        w_train, w_val, w_test
    )
