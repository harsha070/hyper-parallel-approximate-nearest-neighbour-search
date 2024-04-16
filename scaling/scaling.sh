#!/bin/bash

#SBATCH --job-name=scaling
#SBATCH --output=scaling_output/scaling_%Y-%m-%d_%H-%M-%S.out
#SBATCH --error=scaling_output/scaling_%Y-%m-%d_%H-%M-%S.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:10:00
#SBATCH --mail-type=END
#SBATCH --mail-user=larsankile@g.harvard.edu


# Load required modules (if any)
module load gcc/12.1.0-fasrc01

echo "Compiling strong scaling"
# Compile the code
make strong

echo "Running strong scaling"
# Run the code
srun ./strong

# Clean up
rm strong

echo "Plotting strong scaling"
# Plot the results
python plot_strong_scaling.py

echo "Compiling weak scaling"
# Compile the code
make weak

echo "Running weak scaling"
# Run the code
srun ./weak

# Clean up
rm weak

echo "Running weak scaling"
# Plot the results
python plot_weak_scaling.py