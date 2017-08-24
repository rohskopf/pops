import xml.etree.ElementTree as ET
from os import listdir
from os.path import isfile, join
from os import walk

class XMLParser():

    # Initializer
    def __init__(self):
        pass

    # print() overloader
    def __repr__(self):
        return("input!")


    def getConfigQuantities(self, tree):
        """ 
        This function gets quantities associated with a single config
        Inputs: tree - A parsed XML element tree
        Outputs: A list of tuples organized like so:
                 configQuantities is a tuple of quantities associated with config i
                 configQuantities[0] - box
                 configQuantities[1] - atom types
                 configQuantities[2] - positions
                 configQuantities[3] - forces
                 configQuantities[4] - energy
        """
        
        root = tree.getroot()
        # Get atom types
        configAtomSymbols = []
        configAtomTypes = []
        atominfo = root.find('atominfo')
        for array in atominfo.iter('array'):
            name = array.attrib['name']
            if name == 'atoms':
                Set = array.find('set')
                for rc in Set.iter('rc'):
                    c_counter = 1
                    for c in rc.iter('c'):
                        if c_counter==1:
                            configAtomSymbols.append(c.text)
                        if c_counter==2: 
                            val = int(c.text)
                            configAtomTypes.append(val)
                        c_counter = c_counter+1


        configPositions = []
        configForces = []
        configBox = []
        for calculation in root.iter('calculation'):
            calc = calculation.tag
            energies = calculation.findall('energy')

            # Find positions and forces
            for varray in calculation.iter('varray'):
                name = varray.attrib['name']
                tag = varray.tag
                # Positions
                if name == 'positions':
                    #print(name)
                    for v in varray.iter('v'):
                        vtext = v.text
                        vals = [float(a) for a in vtext.split()]
                        configPositions.append(vals)
                # Forces
                if name == 'forces':
                    #print(name)
                    for v in varray.iter('v'):
                        vtext = v.text
                        vals = [float(a) for a in vtext.split()]
                        configForces.append(vals)

            
            # Find final energies
            for energy in energies:
                for i in energy.iter('i'):
                    name = i.attrib['name']
                    if name == 'e_wo_entrp':
                        itext = i.text
                        vals = [float(a) for a in itext.split()]
                        configEnergy = vals[0]

            # Box
            for structure in calculation.iter('structure'):
                for crystal in structure.iter('crystal'):
                    for varray in crystal.iter('varray'):
                        name = varray.attrib['name']
                        if name == 'basis':
                            for v in varray.iter('v'):
                                vtext = v.text
                                vals = [float(a) for a in vtext.split()]
                                configBox.append(vals)

            return (configBox, configAtomTypes, configPositions, configForces, configEnergy)


    def getAllQuantities(self):
        """ 
        This function gets quantities associated with all configs
        Inputs: none
        Outputs: A list of tuples organized like so:
                 allQuantities[i] is a list of quantity i for all configs
                 allQuantities[i][k] is quantity i for config k
                 allQuantities[0] - all boxes
                 allQuantities[1] - all atom types
                 allQuantities[2] - all positions
                 allQuantities[3] - all forces
                 allQuantities[4] - all energies
        """

        # Find all XML files in the /data directory
        files = []
        for (dirpath, dirnames, filenames) in walk('data'):
            if filenames != []:
                files.append(dirpath + '/' + filenames[0])

        numFiles = len(files)
        print('Found %d XML files' % (numFiles))
        # Loop through XML files and get config quantities
        allQuantities = []
        allBox = []
        allAtomTypes = []
        allPositions = []
        allForces = []
        allEnergies = []
        for f in range(0,numFiles):
            #if f >= 1:
            #    break
            if f % 100 == 0:
                print('Parsing %dth XML file' % (f))
            #print(files[f])
            tree = ET.parse(files[f])
            # Get quantities associated with the config here
            configQuantities = self.getConfigQuantities(tree)
            allBox.append(configQuantities[0])
            allAtomTypes.append(configQuantities[1])
            allPositions.append(configQuantities[2])
            allForces.append(configQuantities[3])
            allEnergies.append(configQuantities[4])

        allQuantities = [allBox, allAtomTypes, allPositions, allForces, allEnergies]
        return allQuantities

