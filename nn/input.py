import numpy as np
from descriptors import Descriptors
from ase import Atoms
from ase import units
#from ase.calculators.nn import NeuralNet
from ase.visualize import view
from ase.neighborlist2 import NeighborList

"""
This class prepares everything based on user input and CONFIGS
"""

class Input():
    # Initializer
    def __init__(self):
        pass

    # print() overloader
    def __repr__(self):
        return("input!")

    def readconfigs(self):
        """ Function read CONFIGS file """
        fh = open('CONFIGS', 'r')
        l = fh.readline()
        # Loop through entire file
        configcount = 0
        allnums = []
        allpos = []
        allbox = []
        allU = []
        while (l != ''):
            if "CONFIG" in l:
                configcount += 1
                nums = []
                pos = []
                box = []
                while "/" not in l:
                    l = fh.readline()
                    # Take the \n out
                    if "/" not in l:
                        l = l[:-2]
                        #print(l)
                        vals = [float(x) for x in l.split()]
                        num = vals[0]
                        #print ("num: %d" % num)
                        nums.append(num)
                        x,y,z = vals[1], vals[2], vals[3]
                        pos.append((x,y,z))
                allnums.append(nums)
                allpos.append(pos)
                # Now get the cell 
                l = fh.readline() # potential energy
                U = float(l)
                allU.append(U)
                l = fh.readline() # stress tensor
                l = fh.readline() # total pressure
                l = fh.readline() # box
                vals = [float(x) for x in l.split()]
                xx,yy,zz,xy,xz,yz = vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]
                box.append((xx,0,0))
                box.append((xy,yy,0))
                box.append((xz,yz,zz))
                allbox.append(box)

            l = fh.readline()

        return (allnums, allpos, allbox, allU)

    def storeQuantities(self):
        stuff = self.readconfigs()
        allU = stuff[3]
        return (allU)

    def storeAtoms(self):
        stuff = self.readconfigs()
        allnums = stuff[0]
        allpos = stuff[1]
        allbox = stuff[2]
        numConfigs = len(allpos)
        print("Found %d configs" % (numConfigs))

        allAtoms = []
        for c in range(0,numConfigs):
            # Set up geometry for cth configuration
            config = c
            numbers = allnums[config]
            pos = allpos[config]
            box = allbox[config]
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
        D = Descriptors(atype=1, etas=etas)
        # Get number of configs
        M = len(allAtoms)
        # Loop through configs
        print("Calculating descriptors for all configs...")
        allGList = []
        for m in range(0,M):
            atoms = allAtoms[m]
            configGList = D.calculate_system(atoms, rc)
            allGList.append(configGList)

        print("Calculated descriptors for all configs")
        return allGList 
