# Make an OUTCAR directory to house all the OUTCAR files
mkdir outcar_house

for i in 4.75 5.0 5.1 5.2 5.3 5.4 5.45 5.5 5.55 5.6 5.65 5.7 5.71 5.72 5.73 5.74 5.75 5.755 5.756 5.757 5.758 5.759 5.760 5.761 5.762 5.763 5.764 5.77 5.78 5.79 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7
do

  # Move into the config directory
  cd config${i}

  # Make a copy of OUTCAR called OUTCAR${i}
  cp OUTCAR OUTCAR${i}

  # Copy OUTCAR${i} to /nv/hp13/arohskopf3/scratch/dft_configs/Si/8atoms/outcar_house as "OUTCAR${i}"
  cp OUTCAR${i} /nv/hp13/arohskopf3/scratch/dft_configs/Ge/8atoms/moduli/outcar_house
  
  # Exit the directory and repeat
  cd ..

done
