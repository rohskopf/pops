# Enter path to the LAMMPS src directory
LAMMPS = $(HOME)/lammps/lammps-17Nov16/src

# Enter path to the shared library (.so) file
LAMMPSLIB = ~/lammps/lammps-17Nov16/src/liblammps_g++_openmpi.so


ga: ga.cpp
	mpicxx -std=c++11 -Wno-write-strings -o ga ga.cpp -I$(LAMMPS) $(LAMMPSLIB)
	mpicxx -std=c++11 -Wno-write-strings -o calcz calcz.cpp -I$(LAMMPS) $(LAMMPSLIB)