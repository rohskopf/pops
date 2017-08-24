"""
This class builds a feature vector for a single atom.
atype - type of the atom under consideration
rijList - neighborlist of rij distances, along with atom types j
rc = cutoff 
etas = list of eta values for single atom
N = number of atoms in configuration
"""

import numpy as np
from ase.neighborlist import NeighborList

# Inputs will need to be rc, atoms, atomnumber
class Descriptors():
    # Initializer
    def __init__(self, atype, etas, rc):
        self.atype = atype
        self.etas = etas
        self.rc = rc

    # print() overloader
    def __repr__(self):
        return "%s" % (self.atype)

    # Get neighborlist
    #def nl_calc(self, rc, atoms, 

    # Cutoff calculator
    def fc_calc(self,rij):
        rc = self.rc
        if (rij <= rc):
            fc = 0.5*(np.cos(np.pi*rij/rc) + 1)
        else:
            fc = 0.
        return fc

    # Cutoff list calculator
    def fcList_calc(self,rijList):
        rc = self.rc
        fcList = []
        for i in range(0, len(rijList)):
            fcList.append(self.fc_calc(rijList[i]))
        return fcList

    # G1 calculator for a single atom and single eta
    def G1_calc(self,rijList, eta, fcList):
        rc = self.rc
        G1 = 0.
        for i in range(0, len(rijList)): # loop over neighbors j
            G1 += np.exp((-eta*rijList[i]**2)/(rc**2))*fcList[i]
        return G1

    # G1 list creator (feature vector for a single atom, for all etas)
    def G1List_calc(self,rijList, etas, fcList):
        G1List = []
        for i in range(0, len(etas)): # loop over neighbors j
            G1 = self.G1_calc(rijList, etas[i], fcList)
            G1List.append(G1)
        return G1List

    # Build the feature vector and return it (for a single atom)
    def calculate(self, rijList):
        #cutoffs = [self.rc]
        fcList = self.fcList_calc(rijList)
        G1List = self.G1List_calc(rijList, self.etas, fcList)
        return G1List

    # Build feature vectors for all atoms in the configuration
    def calculate_system(self, atoms):
        rc = self.rc/2. # Need to divide by 2 to get proper cutoff in ASE
        pos = atoms.positions
        N = len(pos)
        cutoffs = [rc] * N # Cutoffs for every atom
        nl = NeighborList(cutoffs=cutoffs, skin=0., sorted=False, self_interaction=False, bothways=True)
        nl.update(atoms)
        configG1List = []
        for a in range(0,N):
            indices, offsets = nl.get_neighbors(a)
            atompos = pos[a]
            #print("atompos: ")
            #print(atompos)
            neighpos = []
            for i, offset in zip(indices, offsets):
                #print(atoms.positions[i] + np.dot(offset, atoms.get_cell()))
                neigh = atoms.positions[i] + np.dot(offset, atoms.get_cell())
                neighpos.append(neigh)
                
                """
                if a == 0:
                #print("neighpos: ")
                #print(neighpos)
                print(len(neighpos))
                """

            """ Now we have the neighbors of the atom, let's loop through calculate descriptors for each neighbor """
            # Calculate fc list and G1_vec
            rijList = []
            #G1 = 0.
            for i in range(0,len(neighpos)):
                """print("Neighbor %d:" % (i) )
                print("neighpos[i][:]: ")
                print(neighpos[i][:])"""
                rij = neighpos[i][:] - atompos
                #print(rij)
                rij = np.linalg.norm(rij)
                #print(rij)
                rijList.append(rij)


            # Make G1List for this group of atoms
            G1List = self.calculate(rijList)
            """
            if a == 0:
                print(neighpos)
                print(rijList)
                print(G1List)
            """
            configG1List.append(G1List)

        # Swap the axes so we have a (MxGxN) list 
        configG1List = np.swapaxes(configG1List, 0, 1)
        #print(np.shape(configG1List))
        return(configG1List)
