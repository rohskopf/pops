module purge
module load intel/15.0
module load fftw/3.3.4
module load mkl/11.2
module load openmpi/1.8
module load octave/3.6.1
module load vasp/5.3.3

# Make sure all the config${i} files are generated
for i in 4.75 5.0 5.1 5.2 5.3 5.4 5.45 5.5 5.55 5.6 5.65 5.7 5.71 5.72 5.73 5.74 5.75 5.755 5.756 5.757 5.758 5.759 5.760 5.761 5.762 5.763 5.764 5.77 5.78 5.79 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7
do

  # Move into the config${i} directory
  cd config${i}

  # Submit a vasp job
  qsub runvasp.pbs

  # Exit the directory and repeat for the next configuration
  cd ..

done
