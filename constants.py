"""
Physical constants and simulation parameters.
"""

# =============================================================================
# Physical Constants
# =============================================================================
SIGMA = 5.67e-8               # Stefan-Boltzmann constant [W/(m²·K⁴)]
BOLTZMANN = 1.380649e-23      # Boltzmann constant [J/K]
ELECTRON_CHARGE = 1.60217663e-19  # Elementary charge [C]


# =============================================================================
# Panel Parameters
# =============================================================================
PANEL_AREA = 1.61             # Panel area [m²]
PCM_THICKNESS = 0.035         # PCM layer thickness [m]


# =============================================================================
# Heat Transfer Coefficients
# =============================================================================
H_CONV_GLASS = 24.486272      # Convection coefficient (glass) [W/(m²·K)]
H_CONV_BACK = 24.486272       # Convection coefficient (layer f) [W/(m²·K)]


# =============================================================================
# PV Panel Specifications (Reference Conditions)
# =============================================================================
G_REF = 1000                  # Reference irradiance [W/m²]
T_REF = 298.13                # Reference temperature [K]
E_G = 1.1                     # Band gap energy [eV]
K_I = 0.005254                # Temperature coefficient of short-circuit current [A/K]

I_SC = 9.06                   # Short-circuit current [A]
V_OC = 46.22                  # Open-circuit voltage [V]
N_S = 72                      # Number of cells in series
N_IDEALITY = 1.3              # Diode ideality factor


# =============================================================================
# Simulation Parameters
# =============================================================================
TIME_STEP = 120               # Time step for temperature derivative [s]
IRRADIANCE_FACTOR = 116       # Irradiance conversion factor

