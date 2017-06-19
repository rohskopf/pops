# Enter path to the LAMMPS src directory
LAMMPS = $(HOME)/lammps/lammps-17Nov16/src

# Enter path to the shared library (.so) file
LAMMPSLIB = ~/lammps/lammps-17Nov16/src/liblammps_g++_openmpi.so


ga: ga_old.cpp
	mpicxx -std=c++11 -Wno-write-strings -o ga_old ga_old.cpp -I$(LAMMPS) $(LAMMPSLIB)
	mpicxx -std=c++11 -Wno-write-strings -o calcz calcz.cpp -I$(LAMMPS) $(LAMMPSLIB)
	mpicxx -std=c++11 -Wno-write-strings -o calcz_ifc calcz_ifc.cpp -I$(LAMMPS) $(LAMMPSLIB)
	mpicxx -std=c++11 -Wno-write-strings -o ga_ifc ga_ifc.cpp -I$(LAMMPS) $(LAMMPSLIB)
