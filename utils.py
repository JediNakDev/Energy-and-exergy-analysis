"""
Utility functions for data processing and analysis.
"""

from typing import List, Union
import numpy as np


# =============================================================================
# Smoothing Functions
# =============================================================================
def smooth(x: Union[List[float], np.ndarray], box_pts: int) -> np.ndarray:
    """
    Apply moving average smoothing filter.
    
    Args:
        x: Input data array
        box_pts: Number of points for the moving average window
    
    Returns:
        Smoothed data array
    """
    box = np.ones(box_pts) / box_pts
    return np.convolve(x, box, mode='same')


def calculate_error(original_data: Union[List[float], np.ndarray], 
                   smoothed_data: Union[List[float], np.ndarray]) -> np.ndarray:
    """
    Calculate the absolute error between original and smoothed data.
    
    Args:
        original_data: Original data array
        smoothed_data: Smoothed data array
    
    Returns:
        Absolute error array
    
    Raises:
        ValueError: If input arrays have different lengths
    """
    if len(original_data) != len(smoothed_data):
        raise ValueError("Input arrays must have the same length.")
    
    original_data = np.array(original_data)
    smoothed_data = np.array(smoothed_data)
    
    return np.abs(original_data - smoothed_data)


# =============================================================================
# Sky Temperature
# =============================================================================
def calculate_sky_temperature(ambient_temps: List[float]) -> List[float]:
    """
    Calculate effective sky temperature from ambient temperature.
    
    Uses the Swinbank correlation: T_sky = 0.0552 * T_amb^1.5
    
    Args:
        ambient_temps: List of ambient temperatures [K]
    
    Returns:
        List of effective sky temperatures [K]
    """
    return [0.0552 * (t ** 1.5) for t in ambient_temps]


# =============================================================================
# Thermal Conductivity
# =============================================================================
def series_thermal_conductance(kperd1: float, kperd2: float) -> float:
    """
    Calculate effective thermal conductance for two layers in series.
    
    Args:
        kperd1: k/d for first layer [W/(m²·K)]
        kperd2: k/d for second layer [W/(m²·K)]
    
    Returns:
        Effective k/d for combined layers [W/(m²·K)]
    """
    return 1 / ((1 / kperd1) + (1 / kperd2))

