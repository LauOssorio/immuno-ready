import joblib
import numpy as np
from tensorflow.keras.models import load_model
import pandas as pd
from src.config import SAVED_MODELS_PATH

def predict_new_samples_cnn_multimodal_classificator(model_path = SAVED_MODELS_PATH + "cnn_multimodal_class.h5", X_new=None, X_cat_new=None, threshold=0.4):
    """
    Predict classification probabilities and binary labels on new data.
    """

    model = load_model(model_path)
    pred_probs = model.predict([X_new, X_cat_new])
    y_pred_prob = pred_probs.flatten()
    y_pred_label = (y_pred_prob > threshold).astype(int)

    return pd.DataFrame({
        'target_prob': y_pred_prob,
        'target_strength': y_pred_label
    })
