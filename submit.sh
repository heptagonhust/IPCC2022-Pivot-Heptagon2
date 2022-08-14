!/bin/bash
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 64
./build/pivot $1