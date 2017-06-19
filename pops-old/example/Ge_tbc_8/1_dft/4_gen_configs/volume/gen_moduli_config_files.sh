# Make sure all the config${i} files that will be generated are deleted before running this in a directory

for i in 4.75 5.0 5.1 5.2 5.3 5.4 5.45 5.5 5.55 5.6 5.65 5.7 5.71 5.72 5.73 5.74 5.75 5.755 5.756 5.757 5.758 5.759 5.760 5.761 5.762 5.763 5.764 5.77 5.78 5.79 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7
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
