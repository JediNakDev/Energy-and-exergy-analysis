import numpy as np
import matplotlib.pyplot as plt

def smooth(x, box_pts):
    box = np.ones(box_pts) / box_pts
    x_smooth = np.convolve(x, box, mode='same')
    return x_smooth

# Example usage:
original_data = [1, 2, 3, 4, 5]
box_pts = 3
smoothed_data = smooth(original_data, box_pts)

# Calculate the error as the absolute difference between original and smoothed data
error = np.abs(original_data - smoothed_data)

# Create a figure and axis
fig, ax = plt.subplots()

# Plot the smoothed data with vertical error bars
ax.errorbar(range(len(smoothed_data)), smoothed_data, yerr=error, fmt='bo-', label='Smoothed Data with Error')

# Add labels and legend
ax.set_xlabel('Data Point')
ax.set_ylabel('Value')
ax.legend()

# Show the plot
plt.show()
