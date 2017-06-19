# Make an LASTGEN directory to house all the LASTGEN files
mkdir history_house

for ((i=1;i<=50;i++))
do

  # Move into the config directory
  cd trial${i}

  # Make a copy of OUTCAR called OUTCAR${i}
  cp HISTORY HISTORY${i}

  # Copy OUTCAR${i} to /nv/hp13/arohskopf3/scratch/dft_configs/Si/8atoms/outcar_house as "OUTCAR${i}"
  cp HISTORY${i} /nv/hp13/arohskopf3/scratch/fitting/Si/8atoms/tvc_yuk/history_house
  
  # Exit the directory and repeat
  cd ..

done
