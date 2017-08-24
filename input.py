import numpy as np
from descriptors import Descriptors
from xmlparser import XMLParser
from ase import Atoms
from ase import units
#from ase.calculators.nn import NeuralNet
from ase.visualize import view

"""
This class prepares everything based on user input and CONFIGS
"""

class Input():
    # Initializer
    def __init__(self):
        parser = XMLParser()
        allQuantities = parser.getAllQuantities()
        self.allBox = allQuantities[0]
        self.allAtomTypes = allQuantities[1]
        self.allPositions = allQuantities[2]
        self.allForces = allQuantities[3]
        self.allEnergies = allQuantities[4]

    # print() overloader
    def __repr__(self):
        return("input!")


    def storeQuantities(self):
        allEnergies = self.allEnergies
        return (allEnergies)

    def storeAtoms(self):
        allnums = self.allAtomTypes
        allpos = self.allPositions
        allbox = self.allBox
        numConfigs = len(allpos)
        print("Found %d configs" % (numConfigs))

        allAtoms = []
        for c in range(0,numConfigs):
            # Set up geometry for cth configuration
            config = c
            numbers = allnums[config]
            pos = allpos[config]
            box = allbox[config]
            pos = np.array(pos)
            box = np.array(box)
            # Convert direct coordinates to cartesian
            pos = np.matmul(pos, box)
            atoms = Atoms(numbers = numbers,
                  positions= pos, 
                  cell=box,
                  pbc = [True, True, True])

            N = len(pos)
            #print('Created %d atoms' % (N) )
            allAtoms.append(atoms)

        return allAtoms

    def storeDescriptors(self, rc, etas):
        allAtoms = self.storeAtoms()
        D = Descriptors(atype=1, etas=etas, rc=rc)
        # Get number of configs
        M = len(allAtoms)
        #print(M)
        # Loop through configs
        print("Calculating descriptors for all configs...")
        allGList = []
        for m in range(0,M):
            if m % 10 == 0:
                print('Calculating descriptors for config %d' % (m))
            atoms = allAtoms[m]
            configGList = D.calculate_system(atoms)
            #print(np.shape(configGList))
            allGList.append(configGList)

        print("Calculated descriptors for all configs")
        return allGList 
