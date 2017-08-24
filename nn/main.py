import sys
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
rc = 2.0
etas = [0.05, 1., 2., 4., 8., 20., 40., 80.]
G = len(etas)

""" Build the network"""
topology = np.array([G,15,1])
activations = "tl"
nn = Network(topology, activations)

""" Calculate all descriptors for all configs"""
inclass = Input()
allGList = inclass.storeDescriptors(rc, etas) #[num configs][descriptor num][atom num]
inputs_list = allGList
#print(inputs_list[0][1][0])

""" Get energies for all configs """
ref_list = inclass.storeQuantities()

M = len(inputs_list)
ref_shape = np.shape(ref_list)
num_ref = ref_shape[0]
if (num_ref != M):
    sys.exit("Number of energies not equal to number of configs")
# Loop through configs to get number of atoms for each config, for normalization purposes
Nk_list = []
for k in range(0,M):
    config = inputs_list[k]
    config_shape = np.shape(config)
    Nk = config_shape[1]
    Nk_list.append(Nk)
Nk_list = np.array(Nk_list)

# Calculate number of weights
hl = np.size(topology)-1
Nw = 0
for k in range(1,hl+1):
    previousN = topology[k-1]
    Nk = topology[k]
    term = previousN*Nk + Nk
    Nw = Nw + term
print('%d weights' % (Nw) )

# Create inputs and outputs placeholder list
inputs = []
outputs = [] # outputs[k] are the energies of each configuration k
for g in range(0,M):
    inputs.append(tf.placeholder(tf.float32, shape=(G,None) ) )
    outputs.append(nn.makeNet(inputs[g]))

# Calculate errors
targets = tf.placeholder(tf.float32, shape=(None) )
diff = (targets - outputs)/Nk_list
diffsq = tf.square(diff)
rmse = tf.sqrt(tf.reduce_mean(diffsq))
loss = rmse

# Other interesting quantities
mae = tf.reduce_mean(diff)
mpe = tf.reduce_mean(diff/targets)*100
mtarget = tf.reduce_mean(targets)

# Define the training step
#train_step = tf.train.GradientDescentOptimizer(1e-3).minimize(loss)
train_step = tf.train.AdamOptimizer(1e-3).minimize(loss)

# Initialize global variables and the session
sess = tf.Session()
init = tf.global_variables_initializer()
sess.run(init)

"""
print("inputs:")
print(sess.run(inputs, feed_dict={i: d for i, d in zip(inputs, configs_list)}))
print("outputs:")
print(sess.run(outputs, feed_dict={i: d for i, d in zip(inputs, configs_list)}))
print("targets:")
print(sess.run(targets, feed_dict={targets: ref_list}))
"""

# Print the trainable variables
"""
variables_names = [v.name for v in tf.trainable_variables()]
values = sess.run(variables_names)
for k, v in zip(variables_names, values):
    print "Variable: ", k
    print "Shape: ", v.shape
    print v
"""

# Prepare feed dictionary for loss calculation
feed_dict={i: d for i, d in zip(inputs, inputs_list)}
feed_dict[targets] = ref_list
for i in range(5000):
    # training
    sess.run(train_step, feed_dict=feed_dict)
    # Let's look at the weights and biases
    if i % 10 == 0:
        # to visualize the result and improvement
        err = sess.run(loss, feed_dict=feed_dict)
        print(i, err)
        # other errors
        #print(sess.run(mae, feed_dict={xs: x_data, ys: y_data}))
        #print(sess.run(mtarget, feed_dict={xs: x_data, ys: y_data}))

        # Print variables if error below tolerance
        if (err < 0.01):
            print(sess.run(mae, feed_dict=feed_dict))
            print(sess.run(mpe, feed_dict=feed_dict))
            variable_names = [v.name for v in tf.trainable_variables()]
            values = np.array( sess.run(variable_names) )
            #write_weights(topology, values)
            break

