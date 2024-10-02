#!/bin/bash

# March 2023

# Automation of creating input files for Koopman's compliancy testing of material

# Run with ./Koopmans_testing_March2023.sh

make_index () {
	local index_name=$1
	shift
	local -a value_array=("$@")
	local i
	# -A means associative array, -g means create a global variable:
	declare -g -A ${index_name}
	for i in "${!value_array[@]}"; do
		eval ${index_name}["${value_array[$i]}"]=$i
	done
}

Category="B_defect CD_defect E_defect"

state1="neutral"
state2="negative1"

e1="electron_density"
e2="ENERGY4Koopmans"

template="input.inp"
pythonfile="xyzlast.py"
inp="input.inp"

log="output.log"

newlogname="KoopmansENERGY.log"

for cat in ${Category} ; do
    X='A'
    Q='PBE-q3'
    types1B=('B1' 'B2')
    make_index types1B_index "${types1B[@]}"
    types2B=('B3' 'B4')
    make_index types2B_index "${types2B[@]}"
    types1C=('C1' 'C2')
    make_index types1C_index "${types1C[@]}"
    types1D=('D1' 'D2')
    make_index types1D_index "${types1D[@]}"
    types2C=('C3' 'C4')
    make_index types2C_index "${types2C[@]}"
    types2D=('D3' 'D4')
    make_index types2D_index "${types2D[@]}"
    types1E=('E1' 'E2')
    make_index types1E_index "${types1E[@]}"
    types2E=('E3' 'E4')
    make_index types2E_index "${types2E[@]}"

    if [ $cat == 'B_defect' ] ; then
        ATOMS=448
        types="B1 B2 B3 B4"
        kinds="dzp-q3\\
     \&END KIND\\
     \&KIND B\\
       BASIS_SET DZVP-MOLOPT-SR-GTH\\
       POTENTIAL GTH-PBE-q5\\
       BASIS_SET AUX_FIT admm-dzp-q5"

    elif [ $cat == 'CD_defect' ] ; then
        ATOMS=448
        types="C1 C1 D1 D1 C3 C4 D3 D4"
        kinds="dzp-q3\\
     \&END KIND\\
     &KIND B\\
       BASIS_SET DZVP-MOLOPT-SR-GTH\\
       POTENTIAL GTH-PBE-q5\\
       BASIS_SET AUX_FIT admm-dzp-q5\\
     \&END KIND\\
     \&KIND F\\
       BASIS_SET DZVP-MOLOPT-SR-GTH\\
       POTENTIAL GTH-PBE-q4\\
       BASIS_SET AUX_FIT admm-dzp-q4"

    elif [ $cat == 'E_defect' ] ; then
        ATOMS=448
        types="E1 E2 E3 E4"
        kinds="dzp-q3 \\
     \&END KIND\\
     \&KIND E\\
       BASIS_SET DZVP-MOLOPT-SR-GTH\\
       POTENTIAL GTH-PBE-q6\\
       BASIS_SET AUX_FIT admm-dzp-q6"

    fi
    for type in ${types} ; do
        dir=$cat/$type
        NAME=calc_${type}

        ## neutral electron_density
        if [[ " ${types1B[*]} " =~ " ${type} " ]]; then
            i=$((${types1B_index[$type]} + 1))
            xyzfile4e1=B1${i}.xyz

        elif [[ " ${types2N[*]} " =~ " ${type} " ]]; then
            i=$((${types2B_index[$type]} + 1))
            xyzfile4e1=B2${i}.xyz

        elif [[ " ${types1C[*]} " =~ " ${type} " ]]; then
            i=$((${types1C_index[$type]} + 1))
            xyzfile4e1=C1${i}.xyz

        elif [[ " ${types1D[*]} " =~ " ${type} " ]]; then
            i=$((${types1D_index[$type]} + 1))
            xyzfile4e1=D1${i}.xyz

        elif [[ " ${types2C[*]} " =~ " ${type} " ]]; then
            i=$((${types2C_index[$type]} + 1))
            xyzfile4e1=C2${i}.xyz

        elif [[ " ${types2D[*]} " =~ " ${type} " ]]; then
            i=$((${types2D_index[$type]} + 1))
            xyzfile4e1=D2${i}.xyz

        elif [[ " ${types1E[*]} " =~ " ${type} " ]]; then
            i=$((${types1E_index[$type]} + 1))
            xyzfile4e1=E1${i}.xyz

        elif [[ " ${types2E[*]} " =~ " ${type} " ]]; then
            i=$((${types2E_index[$type]} + 1))
            xyzfile4e1=E2${i}.xyz

        fi
        work_dir1=$dir/$state1/$e1
        if [ ! -d $work_dir1 ] ; then
            mkdir $work_dir1
        else
            rm -r $work_dir1/*
        fi
        state_dir1=$dir/$state1
        cd $dir
        cp $xyzfile4e1 $state1/$e1/$xyzfile4e1
        cd ../../
        cp $subfile $work_dir1/$subfile
        sed -e "s/LT_name/${NAME}/g" -e "s/#//g" -e "s/LT_X/${X}/g" -e "s/PBE-LT_Q/${Q}/g" -e "s/dzp-LT_Q/${kinds}/g" -e "s/KT_geoxyz/${xyzfile4e1}/g" -e "s/TL_atoms/${ATOMS}/g" \
                   $template > $work_dir1/$inp
        cd $work_dir1
        sbatch $subfile
        cd ../../../../

        ## negative1 koopmans
        xyzfile4e2=${NAME}-pos-L.xyz
        check=$dir/$state2/$xyzfile4e2
        state_dir2=$dir/$state2
        if [ ! -f $check ] ; then
            pos1=${NAME}-pos-1.xyz
            sed -e "s/TL_atoms/${ATOMS}/g" \
                $pythonfile > $state_dir2/$pythonfile
            cd $state_dir2
            python $pythonfile $pos1
            cd ../../../
        fi
        work_dir2=$dir/$state2/$e2
        if [ ! -d $work_dir2 ] ; then
            mkdir $work_dir2
        else
            rm -r $work_dir2/*
        fi
        sed -e "s/${log}/${newlogname}/g" \
                   $subfile > $work_dir2/$subfile
        cd $state_dir2
        cp $xyzfile4e2 $e2/$xyzfile4e2
        cd ../../../
        sed -e "s/LT_name/${NAME}/g" -e "s/LT_X/${X}/g" -e "s/PBE-LT_Q/${Q}/g" -e "s/dzp-LT_Q/${kinds}/g" -e "s/KT_geoxyz/${xyzfile4e2}/g" -e "s/TL_atoms/${ATOMS}/g" \
                   $template > $work_dir2/$inp
        cd $work_dir2
        sbatch $subfile
        cd ../../../../
    done
done