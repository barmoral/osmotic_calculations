#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --job-name=HP_r1
#SBATCH --output=slurm_codes/HP_r1.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export CONCENTRATION=35m
export POSTFIX=r1
export FORCEK=0.68095403

#shouldn't need to change
export CASE=${RESTRAINT}_${CONCENTRATION}
export PDB_INFILE=${CONCENTRATION}_${POSTFIX}.pdb

python build_sim_params.py 
python osmotic_sim_dispatch.py -n ${CASE} -p ${POSTFIX} -pdb ./structures/$PDB_INFILE -omm ${CASE}_${POSTFIX}_sims -spp sim_param_sets/equil_sim.json sim_param_sets/prod_sim.json -r ${RESTRAINT} -k ${FORCEK} -du nanometer
echo "$CASE $POSTFIX run"

sacct --format=jobid,jobname,cputime,elapsed