#!/bin/bash

# Feb 2024

File="A.log"
subfile="submission.sh"
Prep="restart_prep.py"
Rstrt_inp="input_restart.py"
Python="checking.py"
restart="A.restart"
intermediate="A.inter.restart"
intermediate_saved="A.inter"
counter=1

if grep -q "interface.inp" $subfile ;
then 
        if grep -q "OPTIMIZATION STEP:      2" $File ;
        then
                cp $restart inter_files/${intermediate_saved}_initial.restart
                python3 $Prep
                cp $subfile submission_inp.sh
                cp submission_restart.sh $subfile
                # sbatch chaining.sh
        else
                cp interface.inp inter_files/interface_initial.inp
                python3 $Rstrt_inp
                # sbatch chaining.sh
        fi 
else
        if grep -q "PROGRAM ENDED AT" $File ;
        then
                mkdir slurm_files
                mv *.slurm slurm_files
        else
                mv $restart $intermediate
                python3 $Python
                saved=${intermediate_saved}_${counter}.restart
                cd inter_files
                until [ ! -f $saved ] ;
                do
                        ((counter++))
                        echo Counter: $counter
                        saved=${intermediate_saved}_${counter}.restart
                done
                cd ../
                cp $intermediate $saved
                mv $saved inter_files
                mv $intermediate $restart
                # sbatch chaining.sh
        fi
fi