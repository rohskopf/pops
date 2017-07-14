.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program POPS.
# Please modify the variables properly.
# We recommend to use gcc >= 4.8 with openMPI >= 1.8
#-----------------------------------------------

# Enter the path to your LAMMPS directory
LAMMPSDIR = ${HOME}/lammps/lammps-31Mar17/src

# Enter path to the shared library (.so) file
LAMMPSLIB = ~/lammps/lammps-31Mar17/src/liblammps_g++_openmpi.so

#-----------------------------------------------

CXX = mpiCC
CXXFLAGS = -std=c++11 -Wno-write-strings

INCLUDE = -I${LAMMPSDIR}

CXXL =
LDFLAGS = 

LAPACK =
LIBS =

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = pops

CXXSRC= pops.cpp main.cpp popsmemory.cpp popsinput.cpp popstimer.cpp objfun.cpp ga.cpp 

OBJS= ${CXXSRC:.cpp=.o}

default: pops

all: ${PROG}

pops: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LIBS} ${LAMMPSLIB}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
