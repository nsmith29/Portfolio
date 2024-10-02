#!/bin/bash

# Feb 2024

# Job summission script for ARCHER2 supercomputer

#################################################
# Load the relevant CP2K module
# Ensure OMP_NUM_THREADS is consistent with cpus-per-task above
# Launch the executable

module load cp2k/cp2k-9.1.0

export OMP_NUM_THREADS=8
export OMP_PLACES=cores

# Ensure the cpus-per-task option is propagated to srun commands
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK


srun --hint=nomultithread --distribution=block:block cp2k.psmp -i A.inp >> A.log


