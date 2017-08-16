import numpy as np
import tensorflow as tf

class Network():

    # Initializer
    def __init__(self, topology, configs, xs):
        self.topology = topology
        self.configs = configs
        self.xs = xs

    # print() overloader
    def __repr__(self):
        return "%s" % (self.topology)

    def add_layer(self,configs, inputs, in_size, out_size, activation_function):

        # Define variables for this layer
        W = tf.Variable(tf.random_normal([out_size, in_size]))
        b = tf.Variable(tf.random_normal([out_size,1]))

        # Loop through all configs and make a stacked output
        outList = []
        for i in range(0,configs):
            x = inputs[i]
            Wx = tf.matmul(W,x)
            Wx_plus_b = Wx + b
            outList.append(Wx_plus_b)
        outStacked = tf.stack(outList)

        if activation_function is None:
            outputs = outStacked
        else:
            outputs = activation_function(outStacked)
        return outputs

    def makeNet(self):

        add_layer = self.add_layer
        topology = self.topology
        configs = self.configs
        xs = self.xs
    
        print('%d configs\n' % (configs) )
        # Add the connections
        layerdict = {} # Think of this as a connection dictionary!!!
        cs = {}
        # Add first layer (inputs)
        layerdict[0] = add_layer(configs, xs, topology[0], topology[1], activation_function=tf.tanh)
        rest = topology[1:]
        # Add hidden layers
        Nh = len(rest)-1
        print('%d hidden layers\n' % (Nh) )
        for l in range(0,Nh-1):
            #print('%d %d' % (rest[l], rest[l+1]) )
            layerdict[l+1] = add_layer(configs, layerdict[l], rest[l], rest[l+1], activation_function=tf.tanh)

        # Add output layer
        layerdict[len(rest)-1] = add_layer(configs, layerdict[len(rest)-2], rest[len(rest)-2], rest[len(rest)-1], activation_function=None)

        lastlayer = tf.squeeze(layerdict[len(rest)-1])

        return lastlayer

    def add_layer_eval(self,weights,biases,configs, inputs, in_size, out_size, activation_function):


        # Define variables for this layer
        W = tf.constant(weights)
        b = tf.constant(biases)

        # Loop through all configs and make a stacked output
        outList = []
        for i in range(0,configs):
            x = inputs[i]
            Wx = tf.matmul(W,x)
            Wx_plus_b = Wx + b
            outList.append(Wx_plus_b)
        outStacked = tf.stack(outList)

        if activation_function is None:
            outputs = outStacked
        else:
            outputs = activation_function(outStacked)
        return outputs

    def makeNetEval(self, params):

        add_layer_eval = self.add_layer_eval
        topology = self.topology
        configs = self.configs
        xs = self.xs
    
        print('%d configs\n' % (configs) )
        # Add the connections
        layerdict = {} # Think of this as a connection dictionary!!!
        # Add first layer (inputs)
        weightIndx = 0
        biasIndx = 2
        cParams = params[weightIndx:biasIndx]
        weights = cParams[0]
        biases = cParams[1]
        #print(cParams)
        layerdict[0] = add_layer_eval(weights,biases,configs, xs, topology[0], topology[1], activation_function=tf.tanh)
        rest = topology[1:]
        # Add hidden layers
        Nh = len(rest)-1
        print('%d hidden layers\n' % (Nh) )
        for l in range(0,Nh-1):
            weightIndx = weightIndx+2
            biasIndx = biasIndx+2
            cParams = params[weightIndx:biasIndx]
            weights = cParams[0]
            biases = cParams[1]
            #print(cParams)
            #print('%d %d' % (rest[l], rest[l+1]) )
            layerdict[l+1] = add_layer_eval(weights,biases,configs, layerdict[l], rest[l], rest[l+1], activation_function=tf.tanh)

        # Add output layer
        weightIndx = weightIndx+2
        biasIndx = biasIndx+2
        cParams = params[weightIndx:biasIndx]
        weights = cParams[0]
        biases = cParams[1]
        #print(cParams)
        layerdict[len(rest)] = add_layer_eval(weights,biases,configs, layerdict[len(rest)-2], rest[len(rest)-2], rest[len(rest)-1], activation_function=None)

        lastlayer = tf.squeeze(layerdict[len(rest)])

        return lastlayer
