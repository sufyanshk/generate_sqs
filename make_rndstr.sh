#Copyright 2020 (C) Sufyan M. Shaikh
#This script will convert the atomsk generated POSCAR file into the rndstr.in file format
#!/bin/bash/
sed '6,7d' POSCAR > rndstr.in
sed '1,2d' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
sed 's/^......//' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==4{print "1 0 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==5{print "0 1 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==6{print "0 0 1"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in