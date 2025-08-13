import os
import joblib
from tensorflow.keras.callbacks import EarlyStopping
from models.cnn_multimodal_classifier import CNN_multimodal_class
from src.data_processing.pipeline_prepare_training_set import separate_train_val
import matplotlib.pyplot as plt
from src.config import SAVED_MODELS_PATH

def train_model_cnn_multimodal_classificator(save_path=SAVED_MODELS_PATH+"cnn_multimodal_class.h5"):

    # Load data splits
    X_train, X_pca_val,  \
    X_cat_train, X_cat_val, \
    y_train, y_val, \
    w_train, w_val = separate_train_val()

    model = CNN_multimodal_class()

    model.compile(
        optimizer='adam',
        loss='binary_crossentropy',
        metrics=['accuracy']
    )

    es = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)

    history = model.fit(
    x=[X_train, X_cat_train],
    y=y_train,
    sample_weight=w_train,
    validation_data=([X_pca_val, X_cat_val], y_val, w_val),
    epochs=60,
    callbacks=[es],
    verbose=1
    )

    # Save the trained model
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    model.save(save_path)

    # Plot loss curves
    plt.figure(figsize=(12, 5))
    plt.plot(history.history['loss'], label='Loss training')
    plt.plot(history.history['val_loss'], label='Loss validation')
    plt.title('Loss Curve')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    train_model_cnn_multimodal_classificator()
