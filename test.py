"""
Test script for visualizing smoothing error.
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

from utils import smooth


# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    # Example data
    original_data = [1, 2, 3, 4, 5]
    box_pts = 3
    smoothed_data = smooth(original_data, box_pts)
    
    # Calculate error
    error = np.abs(np.array(original_data) - smoothed_data)
    
    # Create plot
    fig, ax = plt.subplots()
    ax.errorbar(
        range(len(smoothed_data)), 
        smoothed_data, 
        yerr=error, 
        fmt='bo-', 
        label='Smoothed Data with Error'
    )
    
    ax.set_xlabel('Data Point')
    ax.set_ylabel('Value')
    ax.legend()
    
    plt.show()
