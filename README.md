# MCSQS from ATAT
This is just a copy of the Alloy Theoretic Automated Toolkit (ATAT) for my reference.

I have added a file named **check\_sqs.sh** which checks whether the *Perfect_match* for the correlation is
found or not.\
If the *Perfect_match* is found in the **bestcorr.out** file, then the code creates a file named **stopsqs**.\
If the *mcsqs* code is running and when it sees the file **stopsqs**, then it cleanly stops its run.

In the **check\_sqs.sh** file you can define how much you want the *mcsqs* code to keep finding out the sqs'es. The time should be mentioned in seconds.

# Disclaimer
This program is not affiliated with ATAT. You are free to modify it, but do so at your own risk.
I have downloaded the ATAT code from [ATAT's official website](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).