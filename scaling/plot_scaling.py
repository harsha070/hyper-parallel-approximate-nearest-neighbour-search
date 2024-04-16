import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV data into pandas DataFrames
strong_data = pd.read_csv("StrongScaling.csv")
weak_data = pd.read_csv("WeakScaling.csv")

# Set up the plot with two subplots side-by-side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

# Plot strong scaling data
ax1.plot(
    strong_data["n_threads"],
    strong_data["time"],
    marker="o",
    linewidth=2,
    color="C0",
    label="Measured Wall Time (s)",
)

ax1.plot(
    strong_data["n_threads"],
    strong_data["time"][0] / strong_data["n_threads"],
    linestyle="--",
    linewidth=2,
    color="C3",
    label="Ideal Scaling",
)

# Customize the strong scaling plot
ax1.set_xlabel("Number of Threads", fontsize=14)
ax1.set_ylabel("Time (s)", fontsize=14)
ax1.set_title("Strong Scaling", fontsize=18)
ax1.tick_params(axis="both", which="major", labelsize=12)
ax1.grid(True)

# Plot weak scaling data
ax2.plot(
    weak_data["n_threads"],
    weak_data["time"],
    marker="o",
    linewidth=2,
    color="C0",
    label="Measured Wall Time (s)",
)

ax2.axhline(
    weak_data["time"][0],
    linestyle="--",
    linewidth=2,
    color="C3",
    label="Ideal Scaling",
)

# Customize the weak scaling plot
ax2.set_xlabel("Number of Threads", fontsize=14)
ax2.set_title("Weak Scaling", fontsize=18)
ax2.tick_params(axis="both", which="major", labelsize=12)
ax2.grid(True)

# Add a legend to the second plot (shared by both subplots)
ax2.legend(fontsize=12)

# Save the plot as an image file (optional)
plt.savefig("scaling_plots.png", bbox_inches="tight", dpi=300)

# Show the plot
plt.show()
