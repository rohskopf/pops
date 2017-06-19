# Make sure all the config${i} files are generated
for ((i=1;i<=100;i++))
do

  # Move into the trial${i} directory
  cd trial${i}

  # Submit a GA run
  qsub runga.pbs

  # Wait 5 seconds for the random number generator seed time
  #sleep 5

  # Exit the directory and repeat for the next trial
  cd ..

done
