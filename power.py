"""
Solar PV power calculation module.

Calculates current-voltage characteristics and maximum power point
for photovoltaic panels based on the single-diode model.
"""

# =============================================================================
# Imports
# =============================================================================
from typing import List
import numpy as np

from constants import (
    BOLTZMANN, ELECTRON_CHARGE,
    G_REF, T_REF, E_G, K_I,
    I_SC, V_OC, N_S, N_IDEALITY
)


# =============================================================================
# Functions
# =============================================================================
def findI(v: List[float], T: List[float], G: List[float]) -> List[List[float]]:
    """
    Calculate current for given voltage, temperature, and irradiance.
    
    Uses the single-diode model for photovoltaic cells.
    
    Args:
        v: Voltage array [V]
        T: Temperature list [K]
        G: Irradiance list [W/m²]
    
    Returns:
        List of current arrays for each temperature/irradiance pair
    """
    i_list = []
    
    for j in range(len(T)):
        t = T[j]
        g = G[j]
        
        # Reverse saturation current at reference
        i_rs = I_SC / (np.exp((ELECTRON_CHARGE * V_OC) / (N_IDEALITY * N_S * BOLTZMANN * t)) - 1)
        
        # Dark saturation current (temperature adjusted)
        i_o = i_rs * ((t / T_REF) ** 3) * np.exp(
            ((ELECTRON_CHARGE * E_G) / (N_IDEALITY * BOLTZMANN)) * ((1 / T_REF) - (1 / t))
        )
        
        # Photogenerated current
        i_ph = (I_SC + K_I * (t - T_REF)) * (g / G_REF)
        
        # Output current (single-diode equation)
        i = i_ph - i_o * (np.exp((ELECTRON_CHARGE * v) / (N_IDEALITY * BOLTZMANN * N_S * t)) - 1)
        i_list.append(i)
    
    return i_list


def findP(i: List[float], v: List[float]) -> float:
    """
    Find maximum power point from I-V curve.
    
    Args:
        i: Current array [A]
        v: Voltage array [V]
    
    Returns:
        Maximum power [W]
    """
    p_max = 0
    for j in range(len(i)):
        p = i[j] * v[j]
        if p > p_max:
            p_max = p
    return p_max


def power(T: List[float], G: List[float]) -> List[float]:
    """
    Calculate maximum power for each temperature/irradiance condition.
    
    Args:
        T: Temperature list [K]
        G: Irradiance list [W/m²]
    
    Returns:
        List of maximum power values [W]
    """
    V = np.arange(-5, V_OC, 0.1)
    I = findI(V, T, G)
    
    return [findP(current, V) for current in I]
