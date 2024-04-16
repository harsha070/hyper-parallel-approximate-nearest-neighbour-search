import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from glob import glob

filename = sorted(glob("strong_output/strong_scaling_*.csv"))[-1]

print(f"Reading {filename}")

# Read the CSV data into a pandas DataFrame
data = pd.read_csv(filename)

# Calculate mean and standard deviation for each number of threads
data = data.groupby("n_threads").agg({"time": ["mean", "std"]})
data.columns = ["time", "std"]

# Set up the plot
plt.figure(figsize=(8, 6))
plt.plot(
    data["time"],
    linestyle="-",
    linewidth=2,
    color="C0",
    label="Strong Scaling",
)

# Add error bars
plt.errorbar(
    data["time"].index,
    data["time"],
    yerr=data["std"],
    linestyle="None",
    marker="o",
    color="C0",
    capsize=5,
    label="Standard Deviation",
)

# Add a line for ideal scaling
plt.plot(
    data["time"].index,
    data["time"].iloc[0] / data["time"].index,
    linestyle="--",
    linewidth=2,
    color="C1",
    label="Ideal Scaling",
)


# Customize the plot
plt.xlabel("Number of Threads", fontsize=14)
plt.ylabel("Time (s)", fontsize=14)
plt.title("Strong Scaling", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)

# Add a legend
plt.legend(fontsize=12)

# Save the plot as an image file (optional)
plt.savefig(
    f"strong_output/strong_scaling_plot_{datetime.now()}.png",
    bbox_inches="tight",
    dpi=300,
)

# Show the plot
plt.show()
