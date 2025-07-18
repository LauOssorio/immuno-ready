import pandas as pd
import joblib

from src.config import PROCESSED_DATA_PATH


from src.data_processing.iedb_data_cleaning import load_clean_iedb
from src.data_processing.normal_data_cleaning import load_clean_normal
from src.data_processing.target_engineering import create_target_features
from src.data_processing.sequence_tokenizer import generate_matrices_for_dataset

def prepare_training_set():
    normals  = load_clean_normal()
    iedb_df = load_clean_iedb()
    # Combine the IEDB and normal dataframes
    full_df= pd.concat([normals, iedb_df], ignore_index=True)
    target_df = create_target_features(full_df)

    # Feature groups
    categorical_features = ['MHC Restriction - Class', 'mhc_status']

    target_encoded = pd.get_dummies(target_df, columns=categorical_features)

    X_pca_aa = generate_matrices_for_dataset(target_encoded)
    X_categorical = target_encoded.iloc[:, -4:]# Assuming the last 4 columns are categorical features
    Y = target_encoded['target_strength']
    sample_weights =target_encoded['averaged_number_positive_subjects_tested']

    joblib.dump(X_pca_aa, PROCESSED_DATA_PATH +'X_pca_aa.joblib', compress=3)
    joblib.dump(X_categorical, PROCESSED_DATA_PATH + 'X_categorical.joblib', compress=3)
    joblib.dump(Y, PROCESSED_DATA_PATH + 'Y.joblib', compress=3)
    joblib.dump(sample_weights, PROCESSED_DATA_PATH + 'sample_weights.joblib', compress=3)



