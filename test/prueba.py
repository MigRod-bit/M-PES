import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.linspace(0, 10, 100)
y = np.sin(x)

# Main Axis (Primary X-axis at the bottom)
fig, ax1 = plt.subplots(figsize=(10, 10))
ax1.plot(x, y, label="Sin(x)")
ax1.set_xlabel("X-axis 1 (Primary)")
ax1.set_xticks([0, 2, 4, 6, 8, 10])
ax1.set_xticklabels(["Zero", "Two", "Four", "Six", "Eight", "Ten"])

# Create 7 additional x-axes using `twiny()` and adjust their positions
axes = [ax1]  # Store all axes for easy management
for i in range(2, 9):  # X-axes 2 to 8
    new_ax = ax1.twiny()  # Create a twin axis
    new_ax.spines["top"].set_position(("axes", 1 + 0.2 * (i - 2)))  # Adjust vertical position
    new_ax.set_xlim(ax1.get_xlim())  # Match the limits with the main axis
    new_ax.set_xticks([0.5 * i for i in range(1, 21, 4)])  # Example ticks
    new_ax.set_xticklabels([f"Label {j}" for j in range(1, len(new_ax.get_xticks()) + 1)])

    axes.append(new_ax)


# Adjust plot to accommodate all the x-axes
fig.subplots_adjust(top=0.9)  # Increase the top margin to prevent overlap
plt.legend()
plt.show()
