module purge
module load intel/15.0
module load fftw/3.3.4
module load mkl/11.2
module load openmpi/1.8
module load octave/3.6.1
module load vasp/5.3.3

# Make sure all the config${i} files are generated
for i in 5.6 5.62 5.64 5.65 5.66 5.67 5.68
do

  # Move into the config${i} directory
  cd config${i}

  # Submit a vasp job
  qsub runvasp.pbs

  # Exit the directory and repeat for the next configuration
  cd ..

done
