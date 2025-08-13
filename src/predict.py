import joblib
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from tensorflow.keras.models import load_model
from sklearn.metrics import (
    precision_score, recall_score, accuracy_score, f1_score,
    confusion_matrix, roc_curve, auc
)

from src.config import SAVED_MODELS_PATH
from src.data_processing.sequence_tokenizer import AA_index_tokenizer
from src.data_processing.cancer_data_cleaning import load_clean_cancer
from src.data_processing.target_engineering import create_target_features


def predict_new_samples_cnn_multimodal_classificator(tokenizer='AA_index_tokenizer', threshold=0.4):
    """
    Predicts immunogenic classification outcomes for new peptide samples using a pretrained
    multimodal CNN model with fixed architecture.

    Returns:
        pd.DataFrame: DataFrame with predicted probabilities, binary predictions, and true labels.
    """
    print("📦 Loading model and test data...")
    model_path = SAVED_MODELS_PATH + "cnn_multimodal_class.h5"
    model = load_model(model_path)

    cancer_df = load_clean_cancer()
    target_cancer = create_target_features(cancer_df)

    print("🔍 Filtering valid peptide sequences...")
    valid_aa_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]{9,20}$')
    target_cancer = target_cancer[target_cancer['Epitope - Name'].apply(lambda x: bool(valid_aa_pattern.match(str(x))))]

    if tokenizer == 'AA_index_tokenizer':
        print("🧬 Tokenizing peptide sequences...")
        X_cancer = AA_index_tokenizer(target_cancer)

    print("📊 Encoding categorical features...")
    categorical_features = ['MHC Restriction - Class', 'mhc_status']
    expected_ohe_columns = [
        'MHC Restriction - Class_Class I',
        'MHC Restriction - Class_Class II',
        'mhc_status_neg',
        'mhc_status_pos'
    ]

    X_cancer_cat = pd.get_dummies(target_cancer, columns=categorical_features)

    drop_cols = [
        'Epitope - Name', 'target_strength',
        '1st in vivo Process - Process Type',
        'Assay - Qualitative Measurement',
        'averaged_number_positive_subjects_tested'
    ]
    X_cancer_cat = X_cancer_cat.drop(columns=[col for col in drop_cols if col in X_cancer_cat.columns], errors='ignore')

    for col in expected_ohe_columns:
        if col not in X_cancer_cat.columns:
            X_cancer_cat[col] = 0
    X_cancer_cat = X_cancer_cat[expected_ohe_columns]

    Y_cancer = target_cancer['target_strength'].astype(int)

    print("📈 Predicting...")
    pred_probs = model.predict([X_cancer, X_cancer_cat], verbose=1)
    y_pred_prob = pred_probs.flatten()
    y_pred_label = (y_pred_prob > threshold).astype(int)

    # Evaluation
    print("✅ Evaluating performance...")
    precision = precision_score(Y_cancer, y_pred_label)
    recall = recall_score(Y_cancer, y_pred_label)
    accuracy = accuracy_score(Y_cancer, y_pred_label)
    f1 = f1_score(Y_cancer, y_pred_label)

    # ROC AUC
    fpr, tpr, thresholds = roc_curve(Y_cancer, y_pred_prob)
    roc_auc = auc(fpr, tpr)

    print()
    print(f"🔍 Precision: {precision:.2f} — How many predicted positives are actually correct.")
    print(f"🎯 Recall: {recall:.2f} — How many actual positives were correctly identified.")
    print(f"✅ Accuracy: {accuracy:.2f} — Overall proportion of correct predictions.")
    print(f"⚖️ F1 Score: {f1:.2f} — Harmonic mean of precision and recall; balances both.")
    print(f"📈 ROC AUC: {roc_auc:.2f} — Probability the model ranks a random positive above a random negative.")

    # Confusion matrix
    cm = confusion_matrix(Y_cancer, y_pred_label)
    plt.figure(figsize=(6, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.show()

    # ROC curve
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC Curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Chance')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - cnn_multimodal_classifier')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return pd.DataFrame({
        'target_prob': y_pred_prob,
        'target_strength': y_pred_label,
        'true_label': Y_cancer.values
    })
