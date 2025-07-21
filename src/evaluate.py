import joblib
import numpy as np
from tensorflow.keras.models import load_model
from src.data_processing.pipeline_prepare_training_set import separate_training_test

from sklearn.metrics import precision_score, recall_score, accuracy_score, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from src.config import SAVED_MODELS_PATH


def evaluate_model_cnn_multimodal_classificator(model_path = SAVED_MODELS_PATH + "cnn_multimodal_class.h5"):

    # Load model and test data
    model = load_model(model_path)
    _, _, X_pca_test, \
    _, _, X_cat_test, \
    _, _, y_test, \
    _, _, _ = separate_training_test()

    pred_probs = model.predict([X_pca_test, X_cat_test])
    y_pred = (pred_probs.flatten() > 0.4).astype(int)

    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    accuracy = accuracy_score(y_test, y_pred)

    print(f"Precision: {precision:.2f}, Recall: {recall:.2f}, Accuracy: {accuracy:.2f}")

    # Confusion matrix plot
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(6,6))
    sns.heatmap(cm, annot=True, fmt="d", cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.show()

if __name__ == "__main__":
    evaluate_model_cnn_multimodal_classificator()
