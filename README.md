# ATAT installation
(ATAT and Atomsk should be downloaded from their respective websites given in the **Dislcaimer** section.)

I suggest to install `ATAT` in the `/bin/` directory so that all the users can use it. If you don't edit the `makefile` of ATAT then only the user who installed it, can use it.  
Execute these commands sequentially (so that all the users will be able to use `ATAT`):
````shell
bash
cd atat  
sed -i '' '1s/\$(HOME)//' makefile #if OSX  
sed -i '1s/\$(HOME)//' makefile #if linux  
make  
make install  
````
# MCSQS from ATAT
I have added a file named `check_sqs.sh` which checks whether the `Perfect_match` for the correlation is found or not. If the `Perfect_match` is found in the `bestcorr.out` file, then the code creates a file named `stopsqs`. If the `mcsqs` code is running and when it sees the file `stopsqs`, then it cleanly stops its run.

In the `check_sqs.sh` file you can define for how many seconds the `mcsqs` code should keep finding out the sqs'es. (I have put the default as 172800 seconds which is equal to 2 days)

# `rndstr.in` from POSCAR generated by `Atomsk`
For generating supercell of specific dimensions, `mcsqs` requires `rndstr.in` file of the required supercell.  
The file `make_rndstr.in` edits the POSCAR file generated from `atomsk` code and creates the `rndstr.in` file.  
This `rndstr.in` can then be given as input to `mcsqs` to generate the SQS of required dimensions.  
(Note: All the atomic sites `rndstr.in` file should be given the relative concentration of the chemical species, according to the alloy chemistry.)

# Disclaimer
The `ATAT` code can be downloaded from ATAT's [official website](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).  
The `Atomsk` code can be downloaded from Atomsk's [official website](https://atomsk.univ-lille.fr/dl.php).  
ATAT and Atomsk are the property of their respective creators.  
No copyright infringement intended.