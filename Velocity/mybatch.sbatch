#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=ss11524@nyu.edu
#SBATCH --output=Su_Velocity_2_5_10_%j.out



module purge
module load matlab/2019b

RUNDIR=$SCRATCH/LOS_ClusterBL/Simulations/
cd $RUNDIR




matlab -nodisplay -nodesktop -r "run Main.m"
