#!/bin/bash
#SBATCH -p long
#SBATCH -a 1-10
#SBATCH --mem=2g
#SBATCH -o results/test_designs3/run.sbatch.out
#SBATCH -e results/test_designs3/run.sbatch.err
source activate 2017_topology_mining
jug execute --jugdir results/test_designs3/jugdir scripts/score_designs.py scripts/test_designs/ results/test_designs3/