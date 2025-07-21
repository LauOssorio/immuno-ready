from src.data_processing.pipeline_prepare_training_set import separate_training_test
import tensorflow as tf
from tensorflow.keras import layers  # Import layers shorthand
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense
from tensorflow.keras import optimizers

