#!/bin/bash

# March 2023

# Adapted and amalgamated version of cutoff_run_March2023.sh and rel_cutoff_run.sh from CP2K website:
# https://manual.cp2k.org/trunk/methods/dft/cutoff.html

# Run with: ./cutoff_run_March2023.sh

LT_type1=LT_value1

cp2k_bin=cp2k.popt
input_file=LT_input
output_file=LT_output
no_proc_per_calc=2
no_proc_to_use=16

counter=1
max_parallel_calcs=$(expr $no_proc_to_use / $no_proc_per_calc)
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry
    cd $work_dir
    if [ -f $output_file ] ; then
        rm $output_file
    fi
    mpirun -np $no_proc_per_calc $cp2k_bin -o $output_file $input_file &
    cd ..
    mod_test=$(echo "$counter % $max_parallel_calcs" | bc)
    if [ $mod_test -eq 0 ] ; then
        wait
    fi
    counter=$(expr $counter + 1)
done
wait
