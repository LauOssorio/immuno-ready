import pandas as pd
import numpy as np
import joblib
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split


from src.config import PROCESSED_DATA_PATH


from src.data_processing.iedb_data_cleaning import load_clean_iedb
from src.data_processing.normal_data_cleaning import load_clean_normal
from src.data_processing.target_engineering import create_target_features
from src.data_processing.sequence_tokenizer import AA_index_tokenizer

def prepare_training_set(tokenizer='AA_index_tokenizer'):
    """
    Prepares and saves the training dataset for immunogenicity prediction.

    This function performs the following steps:
    1. Loads and cleans positive (IEDB) and negative/control (HLA Ligand Atlas) peptide datasets.
    2. Merges datasets into one DataFrame.
    3. Generates target labels using `create_target_features`.
    4. One-hot encodes categorical variables: 'MHC Restriction - Class' and 'mhc_status'.
    5. Tokenizes peptide sequences using the specified tokenizer (default: AA_index).
    6. Extracts one-hot encoded features from the DataFrame.
    7. Extracts the target label (`target_strength`) and sample weights.
    8. Scales the sample weights to a fixed range using MinMaxScaler (e.g., [0.1, 0.2]).
    9. Saves all preprocessed components as `.joblib` files for reuse.

    Parameters
    ----------
    tokenizer : str, optional
        Tokenizer to use for sequence embedding. Currently supports:
        - 'AA_index_tokenizer' (default)

    Saved Files
    -----------
    - X_pca_aa_<tokenizer>.joblib : PCA-transformed amino acid embeddings
    - X_categorical.joblib : One-hot encoded categorical features
    - Y.joblib : Target labels (binary)
    - scaled_sample_weights.joblib : Normalized sample weights

    Returns
    -------
    None
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

    joblib.dump(X_categorical, PROCESSED_DATA_PATH + 'X_categorical.joblib', compress=3)
    joblib.dump(Y, PROCESSED_DATA_PATH + 'Y.joblib', compress=3)
    joblib.dump(scaled_sample_weights, PROCESSED_DATA_PATH + 'scaled_sample_weights.joblib', compress=3)

    #TOKENIZER CHOICE:

    if tokenizer=='AA_index_tokenizer':
        # Generate amino acid PCA embeddings for each peptide sequence (as a list of matrices)
        X_pca_aa_index = AA_index_tokenizer(target_encoded)
        joblib.dump(X_pca_aa_index, PROCESSED_DATA_PATH + 'X_pca_aa_index_tokenized.joblib', compress=3)







def load_or_create_training_data(tokenizer='AA_index_tokenizer'):
    """
    Loads preprocessed training data from disk if available; otherwise, generates it.

    This function attempts to load four joblib files:
    - X_pca_aa_<tokenizer>.joblib : Tokenized and embedded peptide features
    - X_categorical.joblib : One-hot encoded categorical variables
    - Y.joblib : Target labels
    - scaled_sample_weights.joblib : Normalized sample weights

    If any file is missing, it automatically calls `prepare_training_set()` to generate them.

    Parameters
    ----------
    tokenizer : str, optional
        Name of the tokenizer used to generate amino acid embeddings. Default is 'AA_index_tokenizer'.
        This affects the filename of the PCA-encoded peptide feature file.

    Returns
    -------
    tuple
        A 4-tuple containing:
        - X : list or np.ndarray
            Peptide sequence features (e.g., list of matrices)
        - X_categorical : pd.DataFrame or np.ndarray
            One-hot encoded categorical features
        - Y : pd.Series or np.ndarray
            Target variable (binary classification)
        - scaled_sample_weights : np.ndarray
            Scaled sample weights to balance contribution during training
    """
    if tokenizer == "AA_index_tokenizer":
        tokenized_data = 'pca_aa_index_tokenized'

    files = {
        "X_tokenized": PROCESSED_DATA_PATH + "X_" + tokenized_data +".joblib",
        "X_categorical": PROCESSED_DATA_PATH + "X_categorical.joblib",
        "Y": PROCESSED_DATA_PATH + "Y.joblib",
        "sample_weights": PROCESSED_DATA_PATH + "scaled_sample_weights.joblib"
    }

    # Check if all files exist
    if all(os.path.exists(path) for path in files.values()):
        print("✅ All joblib files found. Loading...")
        X = joblib.load(files["X_tokenized"])
        X_categorical = joblib.load(files["X_categorical"])
        Y = joblib.load(files["Y"])
        sample_weights = joblib.load(files["sample_weights"])
    else:
        print("⚠️ One or more files missing. Running prepare_training_set()...")
        prepare_training_set()
        # Now load them
        X = joblib.load(files["X_tokenized"])
        X_categorical = joblib.load(files["X_categorical"])
        Y = joblib.load(files["Y"])
        sample_weights = joblib.load(files["sample_weights"])


    return X, X_categorical, Y, sample_weights




def separate_training_test():
    """
    Loads or generates the processed dataset and splits it into training, validation, and test sets.

    The split is performed in two stages:
    1. 60% of the data is allocated to training.
    2. The remaining 40% is split equally into validation (20%) and test (20%) sets.

    All inputs (peptide embeddings, categorical features, target labels, and sample weights)
    are split in a stratified and consistent way using sklearn's `train_test_split`.

    Returns
    -------
    tuple
        A 12-tuple containing:
        - X_train : list or np.ndarray
            Tokenized and embedded peptide features for training
        - X_val : list or np.ndarray
            Peptide features for validation
        - X_test : list or np.ndarray
            Peptide features for testing
        - X_cat_train : np.ndarray or pd.DataFrame
            One-hot encoded categorical features for training
        - X_cat_val : np.ndarray or pd.DataFrame
            Categorical features for validation
        - X_cat_test : np.ndarray or pd.DataFrame
            Categorical features for testing
        - y_train : np.ndarray or pd.Series
            Target labels for training
        - y_val : np.ndarray or pd.Series
            Target labels for validation
        - y_test : np.ndarray or pd.Series
            Target labels for testing
        - w_train : np.ndarray
            Scaled sample weights for training
        - w_val : np.ndarray
            Sample weights for validation
        - w_test : np.ndarray
            Sample weights for testing
    """
    X, X_categorical, Y, scaled_sample_weights = load_or_create_training_data()

    # First split: Train vs Temp (temp will be split into val and test)
    X_train, X_temp, X_cat_train, X_cat_temp, y_train, y_temp, w_train, w_temp = train_test_split(
        X, X_categorical, Y, scaled_sample_weights,
        test_size=0.4, random_state=42
    )

    # Second split: Validation and Test (from temp)
    X_val, X_test, X_cat_val, X_cat_test, y_val, y_test, w_val, w_test = train_test_split(
        X_temp, X_cat_temp, y_temp, w_temp,
        test_size=0.5, random_state=42
    )

    return (
        X_train, X_val, X_test,
        X_cat_train, X_cat_val, X_cat_test,
        y_train, y_val, y_test,
        w_train, w_val, w_test
    )
