# Make an OUTCAR directory for each c/a ratio to house all the OUTCAR files
mkdir outcar_house

for i in 5.6 5.62 5.64 5.65 5.66 5.67 5.68
do

  # Move into the config directory
  cd config${i}

  # Make a copy of OUTCAR called OUTCAR${i}
  cp OUTCAR OUTCAR${i}

  # Copy OUTCAR${i} to ~/outcar_house as "OUTCAR${i}"
  cp OUTCAR${i} /nv/hp13/arohskopf3/scratch/dft_configs/Ge/8atoms/outcar_house
  
  # Exit the directory and repeat
  cd ..

done
