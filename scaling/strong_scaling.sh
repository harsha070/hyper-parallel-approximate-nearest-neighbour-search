#!/bin/bash

#SBATCH --job-name=strong_scaling
#SBATCH --output=strong_output/strong_scaling_%j.out
#SBATCH --error=strong_output/strong_scaling_%j.err
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
make strong

# Run the code
srun ./strong 32

# Clean up
rm strong

# Plot the results
python plot_strong_scaling.py
