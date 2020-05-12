#Copyright 2020 (C) Sufyan M. Shaikh
#This script will convert the atomsk generated POSCAR file into the rndstr.in file format
#!/bin/bash/

#Be default atomsk creates cartesian coordinates
#But mcsqs requires the rndstr.in file to be in fractional coordinates
#Hence, the below command will convert the cartesian POSCAR into fractional POSCAR
atomsk POSCAR -frac POSCAR

#This will create the rndstr.in file in the required format
#which can be given as input to the mcsqs code
sed '6,7d' POSCAR > rndstr.in
sed '1,2d' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
sed 's/^......//' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==4{print "1 0 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==5{print "0 1 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==6{print "0 0 1"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in