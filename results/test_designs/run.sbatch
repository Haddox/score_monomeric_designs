#!/bin/bash
#SBATCH -p long
#SBATCH -a 1-5
#SBATCH --mem=4g
#SBATCH -o results/test_designs/run.sbatch.out
#SBATCH -e results/test_designs/run.sbatch.err
source activate 2017_topology_mining
jug execute --jugdir results/test_designs/jugdir scripts/score_designs.py scripts/test_designs/ results/test_designs/