#Copyright (C) 2020 S#Copyright (C) 2020 Sufyan M. Shaikh
#!/bin/bash
#PBS -q batch
#PBS -N tihfnbzr_sqs
#PBS -l nodes=node4:ppn=1
#PBS -o tihfnbzr_sqs.out
#PBS -e tihfnbzr_sqs.err

cd $PBS_O_WORKDIR

mcsqs -2=
mcsqs -rc