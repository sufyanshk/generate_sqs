# MCSQS from ATAT
This is just a copy of the Alloy Theoretic Automated Toolkit (ATAT) for my own reference.

I have added a file named **check\_sqs.sh** which checks whether the *Perfect_match* for the correlation is
found or not.\
If the *Perfect_match* is found in the **bestcorr.out** file, then the code creates a file named **stopsqs**.\
If the *mcsqs* code is running and when it sees the file **stopsqs**, then it cleanly stops its run.

In the **check\_sqs.sh** file you can define how much you want the *mcsqs* code to keep finding out the sqs'es. The time should be mentioned in seconds.

### ATAT installation
I suggest to install ATAT in the ````/bin/```` directory so that all the users can use it.\
If you don't edit the **makefile** of ATAT then only the user who installed it, can use it.
````bash
bash
cd atat
(if OSX) sed -i '' '1s/$(HOME)//' makefile
(if linux) sed -i '1s/$(HOME)//' makefile
make
make install
````

# Disclaimer
I have downloaded the ATAT code from [ATAT's official website](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).\
No copyright infringement intended.