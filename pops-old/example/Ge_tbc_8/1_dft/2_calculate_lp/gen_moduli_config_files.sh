# Make sure all the config${i} files that will be generated are deleted before running this in a directory

for i in 5.6 5.62 5.64 5.65 5.66 5.67 5.68
do


  # copy all the files into config${i}
  cp -r files config${i}

  # Move into the config${i} directory
  cd config${i}

  # Change the lattice constant in POSCAR
  cat >POSCAR<<!
  diamond gE:
  $i
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
    8
  direct
  0.25 0.25 0.25
  0.75 0.75 0.25
  0.75 0.25 0.75
  0.25 0.75 0.75
  0.00 0.00 0.00
  0.50 0.50 0.00
  0.50 0.00 0.50
  0.00 0.50 0.50
!
  # Exit the directory and repeat for the next configuration
  cd ..

done
