#!/bin/bash

# March 2023

# Automation of editting CP2K input and cutoff analysis job files for specific calculation materials CP2K cutoff convergence testing is being run for.

# Run with ./Convergence_test_setup_March2023.sh

materialsA="A B"
materialsB="C D E F G H I J K L M N O P"
materialsC="Q R S T U V W X Y Z "


letters="A B C"

cutoffs_var="150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900"
rel_cutoff_var="10 20 30 40 50 60 70 80 90"
rel_cutoff_cst="60"

indicator="rel"

template_file="input.inp"
run_file="cutoff_run.sh"
analysis_file="cutoff_analysis.sh"

# for each group of materials; enter directory of material group
for letter in $letters ; do
    materials = materials{$letter}
    cd $letter

    # for material within group of materials; create directory for material if not already in existence
    for m in $materials ; do
        dir=${m}
        if [ ! -d $dir ] ; then
            mkdir $dir
        fi

        # assign name of certain file-related variables; move xyz file into directory for material
        input_file=${m}.inp
        output_file=${m}.log
        xyz_file=${m}.xyz
        mv $xyz_file $dir/

        # assign variables specifc to material
        if [ ${m} == 'A' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=600
            a=" 5.6456  0.0000  0.0000"
            b=" 0.0000 11.7515  0.0000"
            c=" 0.0000  0.0000  4.3077"
            X=Ab
            Q=q5

        elif [ ${m} == 'B' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=800
            a="10.8690  0.0000  0.0000"
            b=" 0.0000  6.7135  0.0000"
            c=" 0.0000  0.0000  5.5851"
            X=Bb
            Q=q6

        elif [ ${m} == 'C' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=650
            a=" 5.6500  0.0000  0.0000"
            b=" 0.0000  5.6500  0.0000"
            c=" 0.0000  0.0000  5.6500"
            X=C
            Q=q5

        elif [ ${m} == 'D' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=750
            a=" 3.1100  0.0000  0.0000"
            b="-1.5550  2.6933  0.0000"
            c=" 0.0000  0.0000  4.9800"
            X=Da
            Q=q3

        elif [ ${m} == 'E' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=750
            a=" 3.1894  0.0000  0.0000"
            b="-1.5947  2.7621  0.0000"
            c=" 0.0000  0.0000  5.1861"
            X=Ea
            Q=q13

        elif [ ${m} == 'F' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=700
            a=" 3.6153  0.0000  0.0000"
            b=" 0.0000  3.6153  0.0000"
            c=" 0.0000  0.0000  3.6153"
            X=Fa
            Q=q3

        elif [ ${m} == 'G' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=600
            a=" 3.5330  0.0000  0.0000"
            b="-1.7665  3.0597  0.0000"
            c=" 0.0000  0.0000  5.6930"
            X=Ga
            Q=q13

        elif [ ${m} == 'H' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=600
            a=" 3.6800  0.0000  0.0000"
            b="-1.8400  3.1870  0.0000"
            c=" 0.0000  0.0000  6.0100"
            X=Ha
            Q=q13

        elif [ ${m} == 'I' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=800
            a=" 4.4400  0.0000  0.0000"
            b=" 0.0000  4.4400  0.0000"
            c=" 0.0000  0.0000  4.4400"
            X=Ia
            Q=q11

        elif [ ${m} == 'J' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=800
            a=" 9.9676  0.0000  0.0000"
            b=" 0.0000  9.9676  0.0000"
            c=" 0.0000  0.0000  9.9676"
            X=Ja
            Q=q10

        elif [ ${m} == 'K' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=750
            a="11.4670  0.0000  0.0000"
            b=" 0.0000 11.4670  0.0000"
            c=" 0.0000  0.0000 11.4670"
            X=Ka
            Q=q10

        elif [ ${m} == 'L' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=850
            a=" 9.7691  0.0000  0.0000"
            b=" 0.0000  9.7691  0.0000"
            c=" 0.0000  0.0000  9.7691"
            X=La
            Q=q12

        elif [ {m} == 'M' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=850
            a=" 3.8190  0.0000  0.0000"
            b=" 0.0000  3.8190  0.0000"
            c=" 0.0000  0.0000  3.8190"
            X=Ma
            Q=q11

        elif [ ${m} == 'N' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=550
            a=" 4.0290  0.0000  0.0000"
            b="-2.0145  3.4892  0.0000"
            c=" 0.0000  0.0000 22.4250"
            X=Na
            Q=q10

        elif [ ${m} == 'O' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=850
            a=" 2.8413  0.0000  0.0000"
            b="-1.4207  2.4606  0.0000"
            c=" 0.0000  0.0000  9.6930"
            X=Oa
            Q=q4

        elif [ ${m} == 'P' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=850
            a=" 6.7146  0.0000  0.0000"
            b=" 0.0000  7.3091  0.0000"
            c=" 0.0000  0.0000  4.5519"
            X=Pb
            Q=q7
        if [ ${m} == 'Q' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=850
            a=" 4.7657  0.0000  0.0000"
            b="-2.3828  4.1272  0.0000"
            c=" 0.0000  0.0000 13.0100 "
            X=Qa
            Q=q3
        elif [ ${m} == 'R' ] ; then
            if [ $indicator == "rel" ]; then
                cutoffs=700
            a=" 5.0463  0.0000  0.0000 "
            b=" 0.0000  8.7020  0.0000 "
            c=" 0.0000  0.0000  9.2833 "
            X=Ra
            Q=q13
        elif [ ${m} == 'S' ] ; then
            cutoffs=600
            a=" 4.3358  0.0000  0.0000"
            b="-2.1679  3.7549  0.0000"
            c=" 0.0000  0.0000  8.3397 "
            X=Sa
            Q=q3
        elif [ ${m} == 'T' ] ; then
            cutoffs=800
            a="10.1170  0.0000  0.0000"
            b=" 0.0000 10.1170  0.0000"
            c=" 0.0000  0.0000 10.1170"
            X=Ta
            Q=q13
        elif [ ${m} == 'U' ] ; then
            cutoffs=900
            a="10.5700  0.0000  0.0000 "
            b=" 0.0000 10.5700  0.0000 "
            c=" 0.0000  0.0000 10.5700 "
            X=Ua
            Q=q13
        elif [ ${m} == 'V' ] ; then
            cutoffs=800
            a=" 9.8440  0.0000  0.0000 "
            b=" 0.0000  9.8440  0.0000 "
            c=" 0.0000  0.0000  9.844 "
            X=Va
            Q=q11
        elif [ ${m} == 'W' ] ; then
            cutoffs=750
            a=" 4.2170  0.0000  0.0000 "
            b=" 0.0000  4.2170  0.0000 "
            c=" 0.0000  0.0000  4.2170 "
            X=Wa
            Q=q10
        elif [ ${m} == 'X' ] ; then
            cutoffs=550
            a=" 4.8071  0.0000  0.0000 "
            b=" 0.0000  4.8071  0.0000 "
            c=" 0.0000  0.0000  4.8071 "
            X=Xa
            Q=q10
        elif [ ${m} == 'Y' ] ; then
            cutoffs=450
            a=" 3.2499  0.0000  0.0000 "
            b="-1.6249  2.8145  0.0000 "
            c=" 0.0000  0.0000  5.2066 "
            X=Ya
            Q=q12
        elif [ ${m} == 'z' ] ; then
            cutoffs=600
            a=" 4.6837  0.0000 0.0 "
            b=" 0.0000  3.4226 0.0 "
            c="-0.8500  0.0 5.057869394691238 "
            X=Za
            Q=q11
        fi

        # ensure variables are assigned correctly based on indictor.
        if  [ $indicator == "rel" ] ; then
            type1="rel_cutoffs"
            values1=$rel_cutoff_var
            value2=$cutoffs
            type2="cutoff"

        else
            type1="cutoffs"
            values1=$cutoffs_var
            value2=$rel_cutoff_cst
            type2="rel_cutoff"

        # exchange LT_ strings with variables in run_file
        sed -e "s/LT_type1/${type1}/g" \
            -e "s/LT_value1/${values1}/g" \
            -e "s/LT_input/${input_file}/g" \
            -e "s/LT_output/${output_file}/g" \
            $run_file > $dir/$run_file

        # exchange LT strings with variables in analysis file
        sed -e "s/LT_type1/${type1}/g" \
            -e "s/LT_value1/${values1}/g" \
            -e "s/LT_input/${input_file}/g" \
            -e "s/LT_output/${output_file}/g" \
            -e "s/LT_name/${m}/g" \
            -e "s/LT_type2/${type2}/g" \
            -e "s/LT_value2/${value2}/g" \
            $analysis_file > $dir/$analysis_file

        # for each value in values1; create directory for each value if not already in existence
        for ii in $values1 ; do
            work_dir=${dir}/${type1}_${ii}Ry
            if [ ! -d $work_dir ] ; then
                mkdir $work_dir
            else
                rm -r $work_dir/*
            fi

            cp $dir/$xyz_file $work_dir

            if [ $indicator == "rel" ] ; then
                cutoffs=$value2
                rels=${ii}
            else
                cutoffs=${ii}
                rels=$value2

            # exchange LT strings in input file
            sed -e "s/LT_name/${m}/g" \
                -e "s/LT_cutoff/${cutoffs}/g" \
                -e "s/LT_rel_cutoff/${rels}/g" \
                -e "s/LT_a/${a}/g" \
                -e "s/LT_b/${b}/g" \
                -e "s/LT_c/${c}/g" \
                -e "s/LT_geoxyz/${xyz_file}/g" \
                -e "s/LT_x/${X}/g" \
                -e "s/LT_q/${Q}/g" \
                $template_file > $work_dir/$input_file
        done
    done
    cd ../
done