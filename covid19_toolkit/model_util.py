import tensorflow as tf
import numpy as np
import datetime
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense, Conv2DTranspose, Reshape,  Activation
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import plot_model
import numpy as np
import json
import os
import pickle


class Encoder_cnn():
    def __init__(self,
                 conv_filters = [32, 64],
                 conv_kernel_size = [5, 5],
                 conv_stride = [1, 1],
                 input_dim = (30018, 5),
                 latent_dim = 4):

        self.conv_filters = conv_filters
        self.conv_kernel_size = conv_kernel_size
        self.conv_stride = conv_stride
#         self.n_layers = len(conv_filters)
#         self.input_dim = input_dim
        self.latent_dim = latent_dim

    def encoder(self, input_shape):
        # Build 2D convolutional layers
        for i in range(len(self.conv_filters)):
            layers_conv = Conv2D(
                filters = self.conv_filters[i],
                kernel_size = self.conv_kernel_size[i],
                strides = self.conv_stride[i],
                padding = 'same',
                activation = 'relu',
                name = 'cnn_encoder' + str(i)
            )
        # Input
        encoder_input = Input(shape = input_shape[1:], name = 'encoder_input')
#         x = tf.cast(encoder_input, tf.float32)
        
        # Reshape
        x = Reshape(input_shape[1:] + [1])(encoder_input)
        x = tf.cast(x, tf.float32)
        
        # Conv layers
        x = layers_conv(x)

        # Shape
        shape_after_conv = x.shape[1:]

        # Flatten and dense layers
        x = Flatten()(x)
        x = Dense(self.latent_dim, name='encoder_output')(x)

        # Output
        encoder_output = x

        # Model
        encoder_cnn = Model(encoder_input, encoder_output)

        return encoder_cnn, shape_after_conv    


class Decoder_cnn():
    def __init__(self,
                 convTranspose_filters = [64, 32, 1],
                 convTranspose_kernel_size = [5, 5, 5],
                 convTranspose_stride = [1, 1, 1],
                 input_dim = None,
                 latent_dim = 4):

        self.convTranspose_filters = convTranspose_filters
        self.convTranspose_kernel_size = convTranspose_kernel_size
        self.convTranspose_stride = convTranspose_stride
        self.input_dim = input_dim
        
        
    def decoder(self, shape_after_conv):
        
        # Input
        decoder_input = Input(shape = self.input_dim, name = 'decoder_input')
        
        # Dense and Reshape layers
        x = Dense(np.prod(shape_after_conv), name = 'Dense')(decoder_input)
        x = Reshape(shape_after_conv, name = 'Reshape')(x)
        
        # 2D deconvolutional layers
        for i in range(len(self.convTranspose_filters)):
            layers_deconv = Conv2DTranspose(
                filters = self.convTranspose_filters[i], 
                kernel_size = self.convTranspose_kernel_size[i],
                strides = self.convTranspose_stride[i],
                padding = 'same',
#                 activation = 'relu',
                name = 'cnn_decoder' + str(i)
            )
            
            # deconvolutional layers
            x = layers_deconv(x)
            
            # Activation
            if i < len(self.convTranspose_filters) - 1:
                x = ReLU()(x)
            else:
                x = Activation(activation = 'sigmoid')(x)
                
        # Reshape
        x = Reshape((x.shape[1], x.shape[2]))(x)
        
        # Output
        decoder_output = x
        
        # Model
        decoder_cnn = Model(decoder_input, decoder_output)
        
        return decoder_cnn
        
        
# Custom non-fixed input dim               
class Autoencoder_cnn3(tf.keras.Model):        
    def __init__(self, 
                 conv_filters = [32, 64],
                 conv_kernel_size = [5, 5],
                 conv_stride = [1, 1],
                 convTranspose_filters = [64, 32, 1],
                 convTranspose_kernel_size = [5, 5, 5],
                 convTranspose_stride = [1, 1, 1],
                 input_dim = (30018, 5),
                 latent_dim = 4
                ):
        super().__init__()
        self.conv_filters = conv_filters
        self.conv_kernel_size = conv_kernel_size
        self.conv_stride = conv_stride
        self.convTranspose_filters = convTranspose_filters
        self.convTranspose_kernel_size = convTranspose_kernel_size
        self.convTranspose_stride = convTranspose_stride
        self.input_dim = input_dim
        self.latent_dim = latent_dim
    
    
    def build(self, input_shape):
        ecnn = Encoder_cnn(latent_dim = self.latent_dim)
        self.encoder, shape_after_conv = ecnn.encoder(input_shape) 
        dcnn = Decoder_cnn(input_dim = self.latent_dim)
        self.decoder = dcnn.decoder(shape_after_conv)

    
    def call(self, inputs):
        # Autoencoder
        x = self.encoder(inputs)
        x = self.decoder(x)
        return x    
                

        
        
        
        

            
        
        
        
        
        
        