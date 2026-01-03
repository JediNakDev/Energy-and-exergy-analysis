"""
Temperature simulation for conventional PV panel (without PCM).

Layer Structure (top to bottom):
    g: Tempered glass
    b: EVA
    c: PV cells
    d: EVA
    e: Tedlar foil
    f: Transparency acrylic glass

Method:
    - dT/dt = (T[i+1]-T[i])/120
    - When A*X = B, X = invA*B
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

from constants import PANEL_AREA, SIGMA, H_CONV_GLASS, H_CONV_BACK, TIME_STEP, IRRADIANCE_FACTOR
from materials import Material, tempered_glass, eva_front, pv_cells, eva_back, tedlar_foil, acrylic_glass
from utils import smooth, series_thermal_conductance, calculate_sky_temperature
from data_smooth import Tas, Gs
from power import power


# =============================================================================
# Material Instances
# =============================================================================
g = tempered_glass
b = eva_front
c = pv_cells
d = eva_back
e = tedlar_foil
f = acrylic_glass


# =============================================================================
# Simulation Parameters
# =============================================================================
A = PANEL_AREA
h_cong = H_CONV_GLASS
h_conf = H_CONV_BACK

# Derived quantities
G = Gs / IRRADIANCE_FACTOR
Pg = 0.925 * 0.9 * G
T_sky = calculate_sky_temperature(Tas)

# Thermal conductivity between layers
kperdgb = series_thermal_conductance(g.kperd, b.kperd)
kperdbc = series_thermal_conductance(b.kperd, c.kperd)
kperdcd = series_thermal_conductance(c.kperd, d.kperd)
kperdde = series_thermal_conductance(d.kperd, e.kperd)
kperdef = series_thermal_conductance(e.kperd, f.kperd)


# =============================================================================
# Temperature Arrays
# =============================================================================
Teg = [Tas[0], Tas[1]]  # Tempered glass exterior
Tig = [Tas[0], Tas[1]]  # Tempered glass interior
Teb = [Tas[0], Tas[1]]  # EVA exterior
Tib = [Tas[0], Tas[1]]  # EVA interior
Tec = [Tas[0], Tas[1]]  # PV cells exterior
Tic = [Tas[0], Tas[1]]  # PV cells interior
Ted = [Tas[0], Tas[1]]  # EVA exterior
Tid = [Tas[0], Tas[1]]  # EVA interior
Tee = [Tas[0], Tas[1]]  # Tedlar foil exterior
Tie = [Tas[0], Tas[1]]  # Tedlar foil interior
Tef = [Tas[0], Tas[1]]  # Acrylic glass exterior
Tif = [Tas[0], Tas[1]]  # Acrylic glass interior


# =============================================================================
# Main Simulation Loop
# =============================================================================
i = 1
while i < 751:
    idx = 2 * i - 1
    
    # -------------------------------------------------------------------------
    # Radiative Heat Transfer Coefficients
    # -------------------------------------------------------------------------
    h_radag = SIGMA * g.emi * (Teg[idx]**2 + T_sky[idx]**2) * (Teg[idx] + T_sky[idx])
    h_radgb = (SIGMA * (Teb[idx]**2 + Tig[idx]**2) * (Teb[idx] + Tig[idx])) / ((1/g.emi) + (1/b.emi) - 1)
    h_radbc = (SIGMA * (Tec[idx]**2 + Tib[idx]**2) * (Tec[idx] + Tib[idx])) / ((1/b.emi) + (1/c.emi) - 1)
    h_radcd = (SIGMA * (Ted[idx]**2 + Tic[idx]**2) * (Ted[idx] + Tic[idx])) / ((1/c.emi) + (1/d.emi) - 1)
    h_radde = (SIGMA * (Tee[idx]**2 + Tid[idx]**2) * (Tee[idx] + Tid[idx])) / ((1/d.emi) + (1/e.emi) - 1)
    h_radef = (SIGMA * (Tef[idx]**2 + Tie[idx]**2) * (Tef[idx] + Tie[idx])) / ((1/e.emi) + (1/f.emi) - 1)
    h_radfa = SIGMA * f.emi * (Tas[idx]**2 + Tif[idx]**2) * (Tas[idx] + Tif[idx])
    
    h_radg = (SIGMA * (Teg[idx]**2 + Tig[idx]**2) * (Teg[idx] + Tig[idx])) / ((1/g.emi) + (1/g.emi) - 1)
    h_radb = (SIGMA * (Teb[idx]**2 + Tib[idx]**2) * (Teb[idx] + Tib[idx])) / ((1/b.emi) + (1/b.emi) - 1)
    h_radc = (SIGMA * (Tec[idx]**2 + Tic[idx]**2) * (Tec[idx] + Tic[idx])) / ((1/c.emi) + (1/c.emi) - 1)
    h_radd = (SIGMA * (Ted[idx]**2 + Tid[idx]**2) * (Ted[idx] + Tid[idx])) / ((1/d.emi) + (1/d.emi) - 1)
    h_rade = (SIGMA * (Tee[idx]**2 + Tie[idx]**2) * (Tee[idx] + Tie[idx])) / ((1/e.emi) + (1/e.emi) - 1)
    h_radf = (SIGMA * (Tef[idx]**2 + Tif[idx]**2) * (Tef[idx] + Tif[idx])) / ((1/f.emi) + (1/f.emi) - 1)
    
    # -------------------------------------------------------------------------
    # Matrix B Construction (RHS)
    # -------------------------------------------------------------------------
    B0_1 = [-(Pg[2*i]/2) - (h_radag * T_sky[2*i]) - (h_cong * Tas[2*i])]
    B0_2 = [-(Pg[2*i+1]/2) - (h_radag * T_sky[2*i+1]) - (h_cong * Tas[2*i+1])]
    B1_1 = [-(Pg[2*i]/2)]
    B1_2 = [-(Pg[2*i+1]/2)]
    B2, B3, B4, B5, B6, B7, B8, B9, B10 = [0], [0], [0], [0], [0], [0], [0], [0], [0]
    B11_1 = [-(h_radfa * Tas[2*i]) - (h_conf * Tas[2*i])]
    B11_2 = [-(h_radfa * Tas[2*i+1]) - (h_conf * Tas[2*i+1])]
    
    matrixB = np.array([
        B0_1, B1_1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11_1,
        B0_2, B1_2, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11_2
    ])
    
    # -------------------------------------------------------------------------
    # Matrix A Construction (Coefficient matrix)
    # -------------------------------------------------------------------------
    A0 = [-g.kperd - h_radag - h_cong - h_radg, g.kperd + h_radg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    A1 = [-kperdgb + 2*g.kperd + h_radg, -2*g.kperd - h_radg - h_radgb, h_radgb + b.kperd, kperdgb - b.kperd, 0, 0, 0, 0, 0, 0, 0, 0]
    A2 = [kperdgb - g.kperd, g.kperd + h_radgb, -b.kperd - h_radgb - b.kperd - h_radb, -kperdgb + b.kperd + b.kperd + h_radb, 0, 0, 0, 0, 0, 0, 0, 0]
    A3 = [0, 0, -kperdbc + 2*b.kperd + h_radb, -2*b.kperd - h_radb - h_radbc, h_radbc + c.kperd, kperdbc - c.kperd, 0, 0, 0, 0, 0, 0]
    A4 = [0, 0, kperdbc - b.kperd, b.kperd + h_radbc, -c.kperd - h_radbc - c.kperd - h_radc, -kperdbc + c.kperd + c.kperd + h_radc, 0, 0, 0, 0, 0, 0]
    A5 = [0, 0, 0, 0, -kperdcd + 2*c.kperd + h_radc, -2*c.kperd - h_radc - h_radcd, h_radcd + d.kperd, kperdcd - d.kperd, 0, 0, 0, 0]
    A6 = [0, 0, 0, 0, kperdcd - c.kperd, c.kperd + h_radcd, -d.kperd - h_radcd - d.kperd - h_radd, -kperdcd + d.kperd + d.kperd + h_radd, 0, 0, 0, 0]
    A7 = [0, 0, 0, 0, 0, 0, -kperdde + 2*d.kperd + h_radd, -2*d.kperd - h_radd - h_radde, h_radde + e.kperd, kperdde - e.kperd, 0, 0]
    A8 = [0, 0, 0, 0, 0, 0, kperdde - d.kperd, d.kperd + h_radde, -e.kperd - h_radde - e.kperd - h_rade, -kperdde + e.kperd + e.kperd + h_rade, 0, 0]
    A9 = [0, 0, 0, 0, 0, 0, 0, 0, -kperdef + 2*e.kperd + h_rade, -2*e.kperd - h_rade - h_radef, h_radef + f.kperd, kperdef - f.kperd]
    A10 = [0, 0, 0, 0, 0, 0, 0, 0, kperdef - e.kperd, e.kperd + h_radef, -f.kperd - h_radef - f.kperd - h_radf, -kperdef + f.kperd + f.kperd + h_radf]
    A11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, f.kperd + h_radf - h_conf, -f.kperd - h_radf - h_radfa]
    
    vectorOf0 = [0] * 12
    rows = [A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11]
    matrixA = np.array([*[np.append(row, vectorOf0) for row in rows],
                       *[np.append(vectorOf0, row) for row in rows]])
    
    # -------------------------------------------------------------------------
    # Add Time Derivative Terms (dT/dt)
    # -------------------------------------------------------------------------
    difference = np.array([
        [(g.M / (2*A))], [(g.M / (2*A))], [(b.M / (2*A))], [(b.M / (2*A))],
        [(c.M / (2*A))], [(c.M / (2*A))], [(d.M / (2*A))], [(d.M / (2*A))],
        [(e.M / (2*A))], [(e.M / (2*A))], [(f.M / (2*A))], [(f.M / (2*A))],
    ])
    
    for j in range(len(difference)):
        matrixA[j][j] -= difference[j][0] / TIME_STEP
        matrixA[j][j+12] += difference[j][0] / TIME_STEP
        matrixA[j+12][j] -= difference[j][0] / TIME_STEP
        matrixA[j+12][j+12] += difference[j][0] / TIME_STEP
    
    # -------------------------------------------------------------------------
    # Solve System: X = A⁻¹ * B
    # -------------------------------------------------------------------------
    T = np.dot(np.linalg.inv(matrixA), matrixB)
    
    # Append results
    Teg.append(T[0][0]); Tig.append(T[1][0]); Teb.append(T[2][0]); Tib.append(T[3][0])
    Tec.append(T[4][0]); Tic.append(T[5][0]); Ted.append(T[6][0]); Tid.append(T[7][0])
    Tee.append(T[8][0]); Tie.append(T[9][0]); Tef.append(T[10][0]); Tif.append(T[11][0])
    
    Teg.append(T[12][0]); Tig.append(T[13][0]); Teb.append(T[14][0]); Tib.append(T[15][0])
    Tec.append(T[16][0]); Tic.append(T[17][0]); Ted.append(T[18][0]); Tid.append(T[19][0])
    Tee.append(T[20][0]); Tie.append(T[21][0]); Tef.append(T[22][0]); Tif.append(T[23][0])
    
    i += 1


# =============================================================================
# Post-Processing
# =============================================================================
# Calculate PV cell temperature (average of exterior and interior)
Tc_NoPCM = [(Tec[i] + Tic[i]) / 2 for i in range(len(Tec))]

# Smooth temperature data
nTc = np.array(Tc_NoPCM)
Tcs_NoPCM = smooth(nTc - 300, 35) + 300

# Calculate power and energy
P_NoPCM = power(Tc_NoPCM, G)
Ps_NoPCM = power(Tcs_NoPCM, G)
E_NoPCM = [np.trapz(P_NoPCM[:i]) for i in range(len(P_NoPCM))]


# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    print(Gs[330:750].tolist())
