# Make sure all the config${i} files that will be generated are deleted before running this in a directory

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
do

  # copy all the files into config${i}
  cp -r files config${i}

  # Move into the config${i} directory
  cd config${i}

  # Randomize the POSCAR file
  octave randomize.m

  # Submit a vasp job
  #qsub runvasp.pbs

  # Exit the directory and repeat for the next configuration
  cd ..

done
