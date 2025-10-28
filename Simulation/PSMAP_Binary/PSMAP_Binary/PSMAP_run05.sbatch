#!/bin/bash

#SBATCH --output=outputs_%x_%A_%a.out
#SBATCH --error=errors_%x_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mzliu@uchicago.edu

# Load the default version of GNU parallel.
module load parallel

module load gcc
module load rstudio
module load jags

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
ulimit -u 10000

# This specifies the options used to run srun. The "-N1 -n1" options are
# used to allocates a single core to each task.
srun="srun --exclusive -N1 -n1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j $SLURM_NTASKS"

# Use R CMD BATCH to run xxx.R.
$parallel "$srun Rscript run_simulations_cluster.R 0.5 100 0.3 0.16 {1} > output.{1}" ::: {1..10}


