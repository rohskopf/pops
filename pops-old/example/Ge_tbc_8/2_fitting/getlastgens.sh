# Make an LASTGEN directory to house all the LASTGEN files
mkdir lastgen_house

for ((i=1;i<=100;i++))
do

  # Move into the config directory
  cd trial${i}

  # Make a copy of OUTCAR called OUTCAR${i}
  cp LASTGEN LASTGEN${i}

  # Copy OUTCAR${i} to /nv/hp13/arohskopf3/scratch/dft_configs/Si/8atoms/outcar_house as "OUTCAR${i}"
  cp LASTGEN${i} /nv/hp13/arohskopf3/scratch/fitting/Ge/tersoff_8/lastgen_house
  
  # Exit the directory and repeat
  cd ..

done
