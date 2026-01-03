"""
Power output graph for PV panel simulations.

Plots power density profiles comparing PV-PCM and conventional PV systems.
"""

# =============================================================================
# Imports
# =============================================================================
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

from powerdata import NoPCM, P1, P2, P3, P4, P5, P6, P7, M12, M13, M24, M26, M35, M36, M37


# =============================================================================
# Configuration
# =============================================================================
FONT_SIZE_LABELS = 36
FONT_SIZE_TICKS = 32
TIME_START = 150
TIME_END = 570


# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    # Generate time axis
    time_points = []
    for i in range(TIME_START, TIME_END):
        time = datetime(2022, 3, 1, i // 30, (i % 30) * 2)
        time_points.append(time)
    
    # Calculate power differences from conventional PV
    dP1 = [P1[i] - NoPCM[i] for i in range(420)]
    dP2 = [P2[i] - NoPCM[i] for i in range(420)]
    dP3 = [P3[i] - NoPCM[i] for i in range(420)]
    dP4 = [P4[i] - NoPCM[i] for i in range(420)]
    dP5 = [P5[i] - NoPCM[i] for i in range(420)]
    dP6 = [P6[i] - NoPCM[i] for i in range(420)]
    dP7 = [P7[i] - NoPCM[i] for i in range(420)]
    dP12 = [M12[i] - NoPCM[i] for i in range(420)]
    
    # Create plot
    fig, ax = plt.subplots()
    
    # Plot power profiles
    ax.plot_date(time_points, P1, 'k', xdate=True)
    ax.plot_date(time_points, NoPCM, 'k-.', xdate=True)
    
    # Alternative plots (commented out for reference)
    # ax.plot_date(time_points, dP1, 'k', xdate=True)
    # ax.plot_date(time_points, dP2, 'k', dashes=[1, 1.5], xdate=True)
    # ax.plot_date(time_points, dP3, 'k', dashes=[4, 1], xdate=True)
    # ax.plot_date(time_points, dP4, 'k-.', xdate=True)
    # ax.plot_date(time_points, dP5, 'gray', xdate=True)
    # ax.plot_date(time_points, dP6, 'gray', dashes=[1, 1.5], xdate=True)
    # ax.plot_date(time_points, dP7, 'gray', dashes=[4, 1], xdate=True)
    
    # Labels and formatting
    ax.set_ylabel('Power area density [W/$m^{2}$]', fontsize=FONT_SIZE_LABELS)
    ax.set_xlabel('Time', fontsize=FONT_SIZE_LABELS)
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICKS)
    
    # Time format
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    # Limits and grid
    plt.ylim(0, 250)
    plt.grid(True, alpha=0.5)
    
    plt.show()
