from tensorflow.keras import Input, Model, layers

def paper_CNN_multimodal_class(image_shape=(25, 20, 1), categorical_input_shape=(4,)):
    # Image input branch
    image_input = Input(shape=image_shape, name="image_input")
    x = layers.Conv2D(filters=16, kernel_size=(2, 12))(image_input)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    x = layers.Dropout(0.4)(x)

    x = layers.Conv2D(filters=32, kernel_size=(2, 1))(x)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(x)

    x = layers.Flatten()(x)
    x = layers.Dense(128, activation='relu')(x)
    x = layers.Dropout(0.4)(x)

    # Categorical input branch
    cat_input = Input(shape=categorical_input_shape, name="categorical_input")
    y = layers.Dense(32, activation='relu')(cat_input)
    y = layers.Dropout(0.3)(y)

    # Combine both branches
    combined = layers.concatenate([x, y])

    # Output for binary classification
    output = layers.Dense(1, activation='sigmoid', name='classification_output')(combined)

    model = Model(inputs=[image_input, cat_input], outputs=output)
    return model
