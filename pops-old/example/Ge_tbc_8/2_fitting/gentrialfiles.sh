# Make sure all the config${i} files that will be generated are deleted before running this in a directory

for ((i=1;i<=100;i++))
do

  # copy all the files into trial${i}
  cp -r files trial${i}

done
