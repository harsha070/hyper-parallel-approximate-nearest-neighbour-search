#!/bin/bash

#SBATCH --job-name=weak_scaling
#SBATCH --output=weak_output/weak_scaling_%j.out
#SBATCH --error=weak_output/weak_scaling_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=00:10:00
#SBATCH --mail-type=END
#SBATCH --mail-user=larsankile@g.harvard.edu

# Load required modules (if any)
module load gcc/12.1.0-fasrc01

# Compile the code
make weak

# Run the code
srun ./weak

# Clean up
rm weak

# Plot the results
python plot_weak_scaling.py
