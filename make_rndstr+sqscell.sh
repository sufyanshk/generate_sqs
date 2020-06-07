#Copyright 2020 (C) Sufyan M. Shaikh
#This script will convert the atomsk generated POSCAR file into the rndstr.in file format
#!/bin/bash/

#Be default atomsk creates cartesian coordinates
#But mcsqs requires the rndstr.in file to be in fractional coordinates
#Hence, the below command will convert the cartesian POSCAR into fractional POSCAR

#First convert POSCAR to POSCAR.xsf (with fractional coordinates) since
#atomsk gives error if the input file is not named POSCAR
atomsk POSCAR -frac xsf
#Delete the earlier POSCAR
rm POSCAR
#Convert POSCAR.xsf to POSCAR (vasp format)
atomsk POSCAR.xsf pos
#Delete the earlier POSCAR.xsf (which is not requried anymore)
rm POSCAR.xsf

#This will create the rndstr.in file in the required format
#which can be given as input to the mcsqs code
#Delete the line 6-7 (the lines containing no. of atoms and "Direct")
sed '6,7d' POSCAR > rndstr.in
#Delete the lines containing the comment and the scaling factor
sed '1,2d' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
#Delte the intial spaces form the rndstr.in file
sed 's/^.....//' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
#Write the usual lattice vector (these are the multiplication factors in x,y,z directions)
#This multiplication is 100,010 and 001 in x,y and z directions becaue we have already
#Given the actual lattice parameters in the first 3 lines of the rndstr.in file
awk 'NR==4{print "1 0 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==5{print "0 1 0"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in
awk 'NR==6{print "0 0 1"}7' rndstr.in > rndstr.in.tmp && mv rndstr.in.tmp rndstr.in

cat > sqscell.out <<!
1

1 0 0
0 1 0
0 0 1
!