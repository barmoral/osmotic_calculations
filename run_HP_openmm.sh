#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=80:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --job-name=HP_original
#SBATCH --output=slurm_codes/HP_original.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env

#Directory where job is started (used later)
export $START_DIR=$PWD
export $PROJECTS_DIR=/projects/bamo6610/ResearchMS/blanca-openmmions/osmotic-calculations/

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export CONCENTRATION=35m
export POSTFIX=original
export FORCEK=0.68095403

#shouldn't need to change
export CASE=${RESTRAINT}_${CONCENTRATION}
export PDB_INFILE=${CONCENTRATION}_${POSTFIX}.pdb

python build_sim_params.py 
python osmotic_sim_dispatch.py -n ${CASE} -p ${POSTFIX} -pdb ./structures/$PDB_INFILE -omm ${CASE}_${POSTFIX}_sims -spp sim_param_sets/equil_sim.json sim_param_sets/prod_sim.json -r ${RESTRAINT} -k ${FORCEK} -du nanometer
echo "$CASE $POSTFIX run"

#copy back output data to starting directory and to projects directory
cp -Rf ./${CASE}_${POSTFIX}_sims $PROJECTS_DIR
cp -Rf ./${CASE}_${POSTFIX} $PROJECTS_DIR

sacct --format=jobid,jobname,cputime,elapsed
