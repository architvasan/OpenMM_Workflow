#!/bin/bash
#PBS -N test_openmm
#PBS -l select=1
#PBS -l walltime=48:00:00
#PBS -q preemptable
#PBS -l filesystems=home:grand:eagle
#PBS -A datascience
#PBS -o logs/
#PBS -e logs/
#PBS -m abe
#PBS -M avasan@anl.gov

module load conda/2023-10-04
conda activate /eagle/datascience/avasan/envs/open_mm

#cd /lus/grand/projects/datascience/avasan/Simulations/8GCY_Apo/run

python equilibrate.py > equilibrate.log
python production.0.py
