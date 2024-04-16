import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from glob import glob

filename = sorted(glob("weak_output/weak_scaling_*.csv"))[-1]

print(f"Reading {filename}")

# Read the CSV data into a pandas DataFrame
data = pd.read_csv(filename)

# Calculate mean and standard deviation for each number of threads
data = data.groupby("n_threads").agg({"time": ["mean", "std"]})
data.columns = ["time", "std"]

# Calculate the weak efficiency
data["efficiency"] = data["time"].iloc[0] / (data["time"] * data.index)

# Set up the plot
plt.figure(figsize=(8, 6))

# Plot the weak efficiency
plt.plot(
    data["efficiency"],
    linestyle="-",
    linewidth=2,
    color="C0",
    label="Weak Scaling",
)

# Add error bars
plt.errorbar(
    data["efficiency"].index,
    data["efficiency"],
    yerr=data["std"],
    linestyle="None",
    marker="o",
    color="C0",
    capsize=5,
    label="Standard Deviation",
)


# Add a horizontal line for ideal scaling efficiency
plt.axhline(
    y=1,
    linestyle="--",
    linewidth=2,
    color="C1",
    label="Ideal Scaling",
)


# Customize the plot
plt.xlabel("Number of Threads", fontsize=14)
plt.ylabel("Efficiency", fontsize=14)
plt.title("Weak Scaling", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)

# Add a legend
plt.legend(fontsize=12)

# Save the plot as an image file (optional)
plt.savefig(
    f"weak_output/weak_scaling_plot_{datetime.now()}.png", bbox_inches="tight", dpi=300
)

# Show the plot
plt.show()
