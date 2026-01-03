"""
Data smoothing and visualization for temperature and light intensity measurements.
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

from data import Ta, G
from utils import smooth, calculate_error


# =============================================================================
# Data Processing Configuration
# =============================================================================
TRIM_START = 36
TRIM_END = 34
SMOOTHING_WINDOW = 100
TEMP_OFFSET = 300
TEMP_SCALE_FACTOR = 1.2


# =============================================================================
# Data Processing
# =============================================================================
# Trim and offset temperature data
Ta_trimmed = Ta[TRIM_START:-TRIM_END]
Ta_offset = [t - TEMP_OFFSET for t in Ta_trimmed]

# Apply smoothing
Ta_smoothed_offset = smooth(Ta_offset, SMOOTHING_WINDOW)
Ta_smoothed = [t * TEMP_SCALE_FACTOR + TEMP_OFFSET for t in Ta_smoothed_offset]
G_smoothed = smooth(G[TRIM_START:-TRIM_END], SMOOTHING_WINDOW)

# Legacy aliases for backward compatibility
Tas = Ta_smoothed
Gs = G_smoothed

# Calculate errors
T_error = calculate_error(Ta_trimmed, Ta_smoothed)
T_error_percent = (T_error / Ta_smoothed) * 100
G_error = calculate_error(G[TRIM_START:-TRIM_END], G_smoothed)
G_error_percent = (G_error / G_smoothed) * 100


# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    # Print error statistics
    print("Temperature error (%):", np.min(T_error_percent), np.max(T_error_percent), np.median(T_error_percent))
    print("Irradiance error (%):", np.min(G_error_percent), np.max(G_error_percent), np.median(G_error_percent))

    # Generate time axis
    time_points = [datetime(2022, 3, 1, i // 30, (i % 30) * 2) for i in range(150, 570)]

    # Plotting configuration
    PLOT_START = 330
    PLOT_END = 750
    FONT_SIZE_TICKS = 20
    FONT_SIZE_LABELS = 36

    # Create figure with dual y-axes
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # Configure tick sizes
    ax1.tick_params(axis='both', labelsize=FONT_SIZE_TICKS)
    ax2.tick_params(axis='y', labelsize=FONT_SIZE_TICKS)

    # Plot light intensity (primary y-axis)
    ax1.plot_date(time_points, G_smoothed[PLOT_START:PLOT_END], 'k', xdate=True)
    ax1.set_ylabel('Light Intensity [lux]', fontsize=FONT_SIZE_LABELS)
    ax1.set_xlabel('Time', fontsize=FONT_SIZE_LABELS)

    # Plot temperature (secondary y-axis)
    ax2.plot_date(time_points, Ta_smoothed[PLOT_START:PLOT_END], 'k:', xdate=True)
    ax2.set_ylabel('Temperature [K]', fontsize=FONT_SIZE_LABELS)

    # Format x-axis as time
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # Add grid
    plt.grid(True, alpha=0.5)

    plt.show()
