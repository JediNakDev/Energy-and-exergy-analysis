"""
I-V and P-V curve visualization for PV panels.

Plots power-voltage characteristics at different temperatures.
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

from power import findI


# =============================================================================
# Configuration
# =============================================================================
IRRADIANCE = 1000  # [W/m²]
VOLTAGE_MIN = -5   # [V]
VOLTAGE_MAX = 55   # [V]
VOLTAGE_STEP = 0.1 # [V]

# Test temperatures [K]
TEMPERATURES = [
    273.15 + 10,   # 10°C
    273.15 + 25,   # 25°C
    273.15 + 40,   # 40°C
    273.15 + 55,   # 55°C
    273.15 + 70,   # 70°C
]

TEMP_LABELS = ['10°C', '25°C', '40°C', '55°C', '70°C']


# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    # Create voltage array
    V = np.arange(VOLTAGE_MIN, VOLTAGE_MAX, VOLTAGE_STEP)
    
    # Irradiance array (same for all temperatures)
    G = [IRRADIANCE] * len(TEMPERATURES)
    
    # Calculate current for each temperature
    I = findI(V, TEMPERATURES, G)
    
    # Calculate power (P = V * I) for each temperature
    P = [V * I[i] for i in range(len(TEMPERATURES))]
    
    # Create plot
    fig, ax = plt.subplots()
    
    for i, power in enumerate(P):
        ax.plot(V, power, label=TEMP_LABELS[i])
    
    # Labels and formatting
    ax.set_xlabel('Voltage (V)')
    ax.set_ylabel('Power (W)')
    ax.set_xlim(VOLTAGE_MIN, VOLTAGE_MAX)
    ax.set_ylim(0, 400)
    ax.legend()
    ax.minorticks_on()
    
    plt.show()
