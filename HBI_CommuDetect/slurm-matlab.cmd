#!/bin/bash
#########################################################################################
##  slurm-matlab.cmd:                                                                   # 
##   Sample SLURM script for running MATLAB job by single node in HPC2021 system        #
##                                                                                      #
#########################################################################################

#SBATCH --job-name=t-matlab
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bianlb@hku.hk


#SBATCH --partition=c_foss_amd2
#SBATCH --qos=normal


#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128         

#SBATCH --mem-per-cpu=4G

#SBATCH --time=2-00:00:00

#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#########################################################################################

cd ${SLURM_SUBMIT_DIR}
OUTFILE=${SLURM_JOB_NAME}.${SLURM_JOBID}

module load matlab

# MATLAB using SLURM assigned all CPUs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "-- DEMO_individual_LBM_real --" >> ${OUTFILE}
matlab -r 'mcc -m DEMO_individual_LBM_real.m; exit'
# mcc -m DEMO_individual_LBM_real.m
time ./DEMO_individual_LBM_real >> ${OUTFILE}

exit 0
# end of example file