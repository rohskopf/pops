import numpy as np
import tensorflow as tf

class Network():

    # Initializer
    def __init__(self, topology, activations):
        self.topology = topology
        self.activations = activations
        # We need to initialize variables here so that we have a fixed amount
        # Create list of variables
        self.w_list = []
        self.b_list = []
        # Loop over connections
        for c in range(0,len(topology)-1):
            in_size = topology[c]
            out_size = topology[c+1]
            weights = tf.Variable(tf.random_normal([out_size, in_size]))
            biases = tf.Variable(tf.random_normal([out_size, 1]))
            self.w_list.append(weights)
            self.b_list.append(biases)
       

    # print() overloader
    def __repr__(self):
        return "%s" % (self.topology)

    # Return energy for an entire config
    def return_energy(self):
        return tf.reduce_sum(self.energies)

    def add_connection(self,weights, biases,inputs, in_size, out_size, activation):

        #print(weights)
        #print(biases)

        Wx = tf.matmul(weights,inputs)
        Wx_plus_b = Wx + biases

        if activation is 'l':
            outputs = Wx_plus_b
        if activation is 't':
            outputs = tf.tanh(Wx_plus_b)
        return outputs

    def makeNet(self, xs):

        topology = self.topology
        activations = self.activations

        add_connection = self.add_connection

        w_list = self.w_list
        b_list = self.b_list
        #print(w_list)
        #print(b_list)
        C = len(topology)-1 # Number of connections
    
        # Make a dictionary of layer output values
        layerdict = {}
        # Add first layer (inputs)
        layerdict[0] = add_connection(w_list[0], b_list[0],xs, topology[0], topology[1], activation=activations[0])
        # Loop through other connections
        for c in range(1,C):
            layerdict[c] = add_connection(w_list[c], b_list[c],layerdict[c-1], topology[c], topology[c+1], activation=activations[c])

        lastlayer = tf.squeeze(layerdict[C-1])
        
        energy = tf.reduce_sum(lastlayer)

        return energy
