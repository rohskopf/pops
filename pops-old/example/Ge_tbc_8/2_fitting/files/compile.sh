module purge
module load gcc/4.9.0
module load openmpi/1.8

LD_LIBRARY_PATH=/nv/hp13/arohskopf3/scratch/fitting/Si/8atoms/tersoff:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/nv/hp13/arohskopf3/scratch/fitting/Si/8atoms/tersoff:$LD_LIBRARY_PATH

mpiCC -std=c++11 -Wno-write-strings ga.cpp -o ga liblammps.so
