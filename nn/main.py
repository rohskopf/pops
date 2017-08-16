import numpy as np
from descriptors import Descriptors
from input import Input
from ase import Atoms
from ase import units
#from ase.calculators.nn import NeuralNet
from ase.visualize import view
from ase.neighborlist2 import NeighborList

import tensorflow as tf
from network import Network

def write_weights(topology, values):

    """
    for a given layer i, values[i] = output weights
                         values[i+1] = output biases
    """
    #print(values)
    # Analyze the topology to get the max number of neurons - build the weight array
    numLayers = topology.shape[0]
    Nw = 0
    Nwpl = np.zeros(numLayers-1, dtype=int)
    for k in range(1, topology.shape[0]):
      num = topology[k-1]*topology[k] + topology[k]
      Nwpl[k-1] = num
      Nw = Nw + num

    weightLayers = topology[:-1]
    maxNeurons = np.amax(weightLayers)
    maxNwpl = np.amax(Nwpl)
    weightArr = np.empty((numLayers-1, maxNwpl,))
    weightArr[:] = np.NAN

    """
    Now we need to loop through weightArr to add weights as follows:
    [outNeuron1_weight, outneuron1_bias, outneuron2_weight, outneuron2_bias, ....
    [outNeuron1_weight, outneuron1_bias, outneuron2_weight, outneuron2_bias, ....
    [...
    where each row is a connection.
    Keep in mind there are numLayers - 1 connections
    """
    # Loop through connections
    numConnections = numLayers - 1
    for i in range(0,numConnections):
        lindx = 2*i
        bindx = lindx+1
        weights = np.transpose(values[lindx])
        biases = values[bindx]
        shape = np.shape(weights)
        numIn = shape[0]
        numOut = shape[1]
        """
        print("Connection %d: %d inputs and %d outputs" % (i+1, numIn, numOut))
        print("weights:")
        print(weights)
        print("biases:")
        print(biases)
        """

        # Loop over input neurons and store output weights
        numWeights = 0
        neuronIndx = 0
        for j in range(0,numIn):
            outWeights = weights[j]
            numOuts = np.size(outWeights)
            #print('numOuts: %d' % (numOuts))
            # Loop through output weights
            for k in range(0,numOuts):
                weightArr[i][k+neuronIndx] = outWeights[k]
                numWeights = numWeights+1
            neuronIndx = neuronIndx+numOut
            #print(neuronIndx)

        biasIndx = numIn*numOut # Where the biases will start
        # Loop through biases
        for j in range(0,numOut):
            weightArr[i][j+biasIndx] = biases[j]
            
            
   
    #print("weightArr:")
    #print(weightArr)

    # Now we can write the weights
    np.savetxt("weights", weightArr, delimiter = " ")

    return None


""" Declare input parameters """
rc = 5.
etas = [0.001, 0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

""" Calculate all descriptors for all configs"""
inclass = Input()
allGList = inclass.storeDescriptors(rc, etas) #[num configs][atom num][descriptor num (length of eta)]
#print(allGList[1][2][3]) # This is the 4th descriptor for the 2nd config, 3rd atom
allGListSwap = np.swapaxes(allGList, 1, 2)
D = allGListSwap

""" Get energies for all configs """
allU = inclass.storeQuantities()

shape = np.shape(D)
M = shape[0]
G = shape[1]
N = shape[2]
print("%d %d %d" % (M,G,N))

""" Declare x and y data """
x_data = D
y_data = allU

# Declare topology
topology = np.array([G,15,1])
# define placeholder for inputs to network
xs = tf.placeholder(tf.float32)
# Build network and get outputs
nn = Network(topology, M, xs)
outputs = nn.makeNet()

# Take the sum of the outputs along the rows
outsum = tf.reduce_sum(outputs, axis=1)

# Calculate number of weights
hl = np.size(topology)-1
Nw = 0
for k in range(1,hl+1):
    previousN = topology[k-1]
    Nk = topology[k]
    term = previousN*Nk + Nk
    Nw = Nw + term
print('%d weights' % (Nw) )


# Define placeholders for target values
ys = tf.placeholder(tf.float32)
ys = tf.squeeze(ys)

# the error between prediction and real data
diff = (ys - outsum)/N
diffsq = tf.square(diff)
rmse = tf.sqrt(tf.reduce_mean(diffsq))
loss = rmse

# Other interesting quantities
mae = tf.reduce_mean(ys - outsum)
maen = tf.reduce_mean(diff)
mtarget = tf.reduce_mean(ys)
mpe = tf.reduce_mean((ys-outsum)/ys)*100

# Define the training step
#train_step = tf.train.GradientDescentOptimizer(1e-3).minimize(loss)
train_step = tf.train.AdamOptimizer(1e-3).minimize(loss)

# Initialize global variables and the session
sess = tf.Session()
init = tf.global_variables_initializer()
sess.run(init)

# Print the trainable variables
"""
variables_names = [v.name for v in tf.trainable_variables()]
values = sess.run(variables_names)
for k, v in zip(variables_names, values):
    print "Variable: ", k
    print "Shape: ", v.shape
    print v
"""

for i in range(5000):
    # training
    sess.run(train_step, feed_dict={xs: x_data, ys: y_data})
    # Let's look at the weights and biases
    if i % 1 == 0:
        # to visualize the result and improvement
        err = sess.run(loss, feed_dict={xs: x_data, ys: y_data})
        print(i, err)

        # Print variables if error below tolerance
        if (err < 0.01):
            # print other stuff
            print(sess.run(mae, feed_dict={xs: x_data, ys: y_data}))
            print(sess.run(mpe, feed_dict={xs: x_data, ys: y_data}))
            variable_names = [v.name for v in tf.trainable_variables()]
            values = np.array( sess.run(variable_names) )
            write_weights(topology, values)
            break
