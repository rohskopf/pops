# Make an OUTCAR directory to house all the OUTCAR files
mkdir outcar_house

for (( i=1; i<=100; i++ ))
do

  # Move into the config directory
  cd config${i}

  # Make a copy of OUTCAR called OUTCAR${i}
  cp OUTCAR OUTCAR${i}

  # Copy OUTCAR${i} to /nv/hp13/arohskopf3/scratch/dft_configs/Si/8atoms/outcar_house as "OUTCAR${i}"
  cp OUTCAR${i} /nv/hp13/arohskopf3/scratch/dft_configs/Ge/8atoms/randodisp_0.15/outcar_house
  
  # Exit the directory and repeat
  cd ..

done
