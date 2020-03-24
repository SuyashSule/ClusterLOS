#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=ss11524@nyu.edu
#SBATCH --output=suyash_%j.out



module purge
module load matlab/2019b

RUNDIR=$SCRATCH/LOS_ClusterBL/Simulations/
cd $RUNDIR




matlab -nodisplay -nodesktop -r "run Main.m"
