rm WAVECAR
for i in 2 4 6 8 11 13 15 17 19
do
cat >KPOINTS <<!
K-Points
 0
Monkhorst Pack
 $i $i $i
 0  0  0
!
echo "a= $i" ; vasp
E=`tail -1 OSZICAR` ; echo $i $E >>SUMMARY.kpoints
octave --silent --eval "getdata('$i')"
done
cat SUMMARY.kpoints
