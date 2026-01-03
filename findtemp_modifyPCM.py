"""
Temperature simulation for PV panel with modified dual-PCM layer.

TODO: Edit PCMs' properties in the PCM configuration section.

Layer Structure (top to bottom):
    g: Tempered glass
    b: EVA
    c: PV cells
    d: EVA
    e: Tedlar foil
    p: PCM (mixture of p1 and p2)
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
import matplotlib.dates as mdates
from datetime import datetime

from constants import PANEL_AREA, SIGMA, H_CONV_GLASS, H_CONV_BACK, TIME_STEP, IRRADIANCE_FACTOR
from materials import (
    Material, PCM, tempered_glass, eva_front, pv_cells, eva_back, 
    tedlar_foil, acrylic_glass, create_pcm
)
from utils import smooth, series_thermal_conductance, calculate_sky_temperature
from data_smooth import Tas, Gs
from power import power
from findtemp_NoPCM import P_NoPCM, Tc_NoPCM, Ps_NoPCM, Tcs_NoPCM


# =============================================================================
# Material Instances
# =============================================================================
g = tempered_glass
b = eva_front
c = pv_cells
d = eva_back
e = tedlar_foil
f = acrylic_glass

# TODO: PCM configuration - change PCM types here
p1 = create_pcm("paraffin_wax")
p2 = create_pcm("sodium_sulfate_decahydrate")


# =============================================================================
# Physical Constants
# =============================================================================
A = PANEL_AREA
h_cong = H_CONV_GLASS
h_conf = H_CONV_BACK


# =============================================================================
# Optimization Loop (PCM ratio sweep)
# =============================================================================
P_avg = []
Cost = []
Cost_eff = []
ratio = range(2, 3)  # PCM1 ratio percentage

for q in ratio:
    r = q / 100  # Convert to fraction
    
    # -------------------------------------------------------------------------
    # Calculate Mixed PCM Properties
    # -------------------------------------------------------------------------
    M1s2s = p1.M_s * r + p2.M_s * (1 - r)   # Both solid
    M1l2s = p1.M_l * r + p2.M_s * (1 - r)   # PCM1 liquid, PCM2 solid
    M1l2l = p1.M_l * r + p2.M_l * (1 - r)   # Both liquid
    M1sl = p1.M_sl * r                       # PCM1 latent heat
    M2sl = p2.M_sl * (1 - r)                 # PCM2 latent heat
    p_emi = p1.emi * r + p2.emi * (1 - r)   # Mixed emissivity
    
    print(M1s2s)
    
    state = [0, 0]  # [PCM1_state, PCM2_state]: 0=solid, 1=liquid
    
    # -------------------------------------------------------------------------
    # Derived Quantities
    # -------------------------------------------------------------------------
    G = Gs / IRRADIANCE_FACTOR
    Pg = 0.925 * 0.9 * G
    Mp = M1s2s
    T_sky = calculate_sky_temperature(Tas)
    
    changestate = False
    
    # Thermal conductivity between layers
    kperdgb = series_thermal_conductance(g.kperd, b.kperd)
    kperdbc = series_thermal_conductance(b.kperd, c.kperd)
    kperdcd = series_thermal_conductance(c.kperd, d.kperd)
    kperdde = series_thermal_conductance(d.kperd, e.kperd)
    
    # -------------------------------------------------------------------------
    # Initialize Temperature Arrays
    # -------------------------------------------------------------------------
    Teg = [Tas[0], Tas[1]]
    Tig = [Tas[0], Tas[1]]
    Teb = [Tas[0], Tas[1]]
    Tib = [Tas[0], Tas[1]]
    Tec = [Tas[0], Tas[1]]
    Tic = [Tas[0], Tas[1]]
    Ted = [Tas[0], Tas[1]]
    Tid = [Tas[0], Tas[1]]
    Tee = [Tas[0], Tas[1]]
    Tie = [Tas[0], Tas[1]]
    Tp = [Tas[0], Tas[1]]
    Tef = [Tas[0], Tas[1]]
    Tif = [Tas[0], Tas[1]]
    
    # -------------------------------------------------------------------------
    # Helper Functions
    # -------------------------------------------------------------------------
    def calculate_radiation_coefficients(idx):
        """Calculate radiative heat transfer coefficients at index."""
        h = {}
        h['h_radag'] = SIGMA * g.emi * (Teg[idx]**2 + T_sky[idx]**2) * (Teg[idx] + T_sky[idx])
        h['h_radgb'] = (SIGMA * (Teb[idx]**2 + Tig[idx]**2) * (Teb[idx] + Tig[idx])) / ((1/g.emi) + (1/b.emi) - 1)
        h['h_radbc'] = (SIGMA * (Tec[idx]**2 + Tib[idx]**2) * (Tec[idx] + Tib[idx])) / ((1/b.emi) + (1/c.emi) - 1)
        h['h_radcd'] = (SIGMA * (Ted[idx]**2 + Tic[idx]**2) * (Ted[idx] + Tic[idx])) / ((1/c.emi) + (1/d.emi) - 1)
        h['h_radde'] = (SIGMA * (Tee[idx]**2 + Tid[idx]**2) * (Tee[idx] + Tid[idx])) / ((1/d.emi) + (1/e.emi) - 1)
        h['h_radep'] = (SIGMA * (Tp[idx]**2 + Tie[idx]**2) * (Tp[idx] + Tie[idx])) / ((1/e.emi) + (1/p_emi) - 1)
        h['h_radpf'] = (SIGMA * (Tef[idx]**2 + Tp[idx]**2) * (Tef[idx] + Tp[idx])) / ((1/p_emi) + (1/f.emi) - 1)
        h['h_radfa'] = SIGMA * f.emi * (Tas[idx]**2 + Tif[idx]**2) * (Tas[idx] + Tif[idx])
        
        h['h_radg'] = (SIGMA * (Teg[idx]**2 + Tig[idx]**2) * (Teg[idx] + Tig[idx])) / ((1/g.emi) + (1/g.emi) - 1)
        h['h_radb'] = (SIGMA * (Teb[idx]**2 + Tib[idx]**2) * (Teb[idx] + Tib[idx])) / ((1/b.emi) + (1/b.emi) - 1)
        h['h_radc'] = (SIGMA * (Tec[idx]**2 + Tic[idx]**2) * (Tec[idx] + Tic[idx])) / ((1/c.emi) + (1/c.emi) - 1)
        h['h_radd'] = (SIGMA * (Ted[idx]**2 + Tid[idx]**2) * (Ted[idx] + Tid[idx])) / ((1/d.emi) + (1/d.emi) - 1)
        h['h_rade'] = (SIGMA * (Tee[idx]**2 + Tie[idx]**2) * (Tee[idx] + Tie[idx])) / ((1/e.emi) + (1/e.emi) - 1)
        h['h_radf'] = (SIGMA * (Tef[idx]**2 + Tif[idx]**2) * (Tef[idx] + Tif[idx])) / ((1/f.emi) + (1/f.emi) - 1)
        return h
    
    def build_matrix_A(h):
        """Build coefficient matrix A."""
        A0 = [-g.kperd - h['h_radag'] - h_cong - h['h_radg'], g.kperd + h['h_radg'], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A1 = [-kperdgb + 2*g.kperd + h['h_radg'], -2*g.kperd - h['h_radg'] - h['h_radgb'], h['h_radgb'] + b.kperd, kperdgb - b.kperd, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A2 = [kperdgb - g.kperd, g.kperd + h['h_radgb'], -b.kperd - h['h_radgb'] - b.kperd - h['h_radb'], -kperdgb + b.kperd + b.kperd + h['h_radb'], 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A3 = [0, 0, -kperdbc + 2*b.kperd + h['h_radb'], -2*b.kperd - h['h_radb'] - h['h_radbc'], h['h_radbc'] + c.kperd, kperdbc - c.kperd, 0, 0, 0, 0, 0, 0, 0]
        A4 = [0, 0, kperdbc - b.kperd, b.kperd + h['h_radbc'], -c.kperd - h['h_radbc'] - c.kperd - h['h_radc'], -kperdbc + c.kperd + c.kperd + h['h_radc'], 0, 0, 0, 0, 0, 0, 0]
        A5 = [0, 0, 0, 0, -kperdcd + 2*c.kperd + h['h_radc'], -2*c.kperd - h['h_radc'] - h['h_radcd'], h['h_radcd'] + d.kperd, kperdcd - d.kperd, 0, 0, 0, 0, 0]
        A6 = [0, 0, 0, 0, kperdcd - c.kperd, c.kperd + h['h_radcd'], -d.kperd - h['h_radcd'] - d.kperd - h['h_radd'], -kperdcd + d.kperd + d.kperd + h['h_radd'], 0, 0, 0, 0, 0]
        A7 = [0, 0, 0, 0, 0, 0, -kperdde + 2*d.kperd + h['h_radd'], -2*d.kperd - h['h_radd'] - h['h_radde'], h['h_radde'] + e.kperd, kperdde - e.kperd, 0, 0, 0]
        A8 = [0, 0, 0, 0, 0, 0, kperdde - d.kperd, d.kperd + h['h_radde'], -e.kperd - h['h_radde'] - e.kperd - h['h_rade'], -kperdde + e.kperd + e.kperd + h['h_rade'], 0, 0, 0]
        A9 = [0, 0, 0, 0, 0, 0, 0, 0, h['h_rade'], -h['h_rade'] - h['h_radep'], h['h_radep'], 0, 0]
        A10 = [0, 0, 0, 0, 0, 0, 0, 0, e.kperd, -e.kperd + h['h_radep'], -h['h_radep'] - h['h_radpf'], f.kperd + h['h_radpf'], -f.kperd]
        A11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, h['h_radpf'], -f.kperd - h['h_radpf'] - f.kperd - h['h_radf'], f.kperd + f.kperd + h['h_radf']]
        A12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, f.kperd + h['h_radf'], -f.kperd - h['h_radf'] - h['h_radfa'] - h_conf]
        
        vectorOf0 = [0] * 13
        rows = [A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12]
        return np.array([*[np.append(row, vectorOf0) for row in rows],
                        *[np.append(vectorOf0, row) for row in rows]])
    
    # -------------------------------------------------------------------------
    # Main Simulation Loop
    # -------------------------------------------------------------------------
    i = 1
    while i < 380:
        
        # Check PCM state
        if Tp[2*i-1] <= p1.Tm and state == [0, 0]:
            changestate = False
            Mp = M1s2s
        elif Tp[2*i-1] >= p2.Tm + 1 and state == [1, 1]:
            changestate = False
            Mp = M1l2l
        elif Tp[2*i-1] >= p1.Tm + 1 and Tp[2*i-1] <= p2.Tm and state == [1, 0]:
            changestate = False
            Mp = M1l2s
        elif i == 1:
            changestate = True
        else:
            if not changestate:
                print(i)
                for arr in [Teg, Tig, Teb, Tib, Tec, Tic, Ted, Tid, Tee, Tie, Tp, Tef, Tif]:
                    arr.pop()
                    arr.pop()
                i -= 1
            changestate = True
        
        # ---------------------------------------------------------------------
        # Normal Temperature Evolution (no phase change)
        # ---------------------------------------------------------------------
        if not changestate:
            idx = 2 * i - 1
            h = calculate_radiation_coefficients(idx)
            matrixA = build_matrix_A(h)
            
            # Build matrix B
            B0_1 = [-(Pg[2*i]/2) - (h['h_radag'] * T_sky[2*i]) - (h_cong * Tas[2*i])]
            B0_2 = [-(Pg[2*i+1]/2) - (h['h_radag'] * T_sky[2*i+1]) - (h_cong * Tas[2*i+1])]
            B1_1, B1_2 = [-(Pg[2*i]/2)], [-(Pg[2*i+1]/2)]
            B2, B3, B4, B5, B6, B7, B8, B9, B10, B11 = [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]
            B12_1 = [-(h['h_radfa'] * Tas[2*i]) - (h_conf * Tas[2*i])]
            B12_2 = [-(h['h_radfa'] * Tas[2*i+1]) - (h_conf * Tas[2*i+1])]
            
            matrixB = np.array([B0_1, B1_1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12_1,
                               B0_2, B1_2, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12_2])
            
            # Add time derivative
            difference = np.array([
                [(g.M / (2*A))], [(g.M / (2*A))], [(b.M / (2*A))], [(b.M / (2*A))],
                [(c.M / (2*A))], [(c.M / (2*A))], [(d.M / (2*A))], [(d.M / (2*A))],
                [(e.M / (2*A))], [(e.M / (2*A))], [(Mp / A)], [(f.M / (2*A))], [(f.M / (2*A))],
            ])
            
            for j in range(len(difference)):
                matrixA[j][j] -= difference[j][0] / TIME_STEP
                matrixA[j][j+13] += difference[j][0] / TIME_STEP
                matrixA[j+13][j] -= difference[j][0] / TIME_STEP
                matrixA[j+13][j+13] += difference[j][0] / TIME_STEP
            
            # Solve system
            T = np.dot(np.linalg.inv(matrixA), matrixB)
            
            # Append results
            for arr, idx_t in [(Teg, 0), (Tig, 1), (Teb, 2), (Tib, 3), (Tec, 4), (Tic, 5),
                              (Ted, 6), (Tid, 7), (Tee, 8), (Tie, 9), (Tp, 10), (Tef, 11), (Tif, 12)]:
                arr.append(T[idx_t][0])
            for arr, idx_t in [(Teg, 13), (Tig, 14), (Teb, 15), (Tib, 16), (Tec, 17), (Tic, 18),
                              (Ted, 19), (Tid, 20), (Tee, 21), (Tie, 22), (Tp, 23), (Tef, 24), (Tif, 25)]:
                arr.append(T[idx_t][0])
            
            i += 1
        
        # ---------------------------------------------------------------------
        # Phase Change State
        # ---------------------------------------------------------------------
        else:
            if (Tp[2*i-1] <= p2.Tm and state == [1, 0]) or state == [0, 0]:
                Mp = M1sl
                Tpchange = p1.Tm + 0.5
                s = 1
            else:
                Mp = M2sl
                Tpchange = p2.Tm + 0.5
                s = 2
            
            totalHeat = 0
            print(i * 2, state, r)
            Tp[(i * 2) - 1] = Tpchange
            
            while totalHeat < Mp and totalHeat >= 0:
                idx = 2 * i - 1
                h = calculate_radiation_coefficients(idx)
                matrixA = build_matrix_A(h)
                
                # Build B matrix for phase change
                print(i, '###')
                B0_1 = [-(Pg[2*i]/2) - (h['h_radag'] * T_sky[2*i]) - (h_cong * Tas[2*i])]
                B0_2 = [-(Pg[2*i+1]/2) - (h['h_radag'] * T_sky[2*i+1]) - (h_cong * Tas[2*i+1])]
                B1_1, B1_2 = [-(Pg[2*i]/2)], [-(Pg[2*i+1]/2)]
                B2, B3, B4, B5, B6, B7, B8 = [0], [0], [0], [0], [0], [0], [0]
                B9 = [-(h['h_radep'] * Tpchange)]
                B11 = [-(h['h_radpf'] * Tpchange)]
                B12_1 = [-(h['h_radfa'] * Tas[2*i]) - (h_conf * Tas[2*i])]
                B12_2 = [-(h['h_radfa'] * Tas[2*i+1]) - (h_conf * Tas[2*i+1])]
                
                matrixB = np.array([B0_1, B1_1, B2, B3, B4, B5, B6, B7, B8, B9, B11, B12_1,
                                   B0_2, B1_2, B2, B3, B4, B5, B6, B7, B8, B9, B11, B12_2])
                
                # Remove PCM rows/columns
                matrixA = np.delete(matrixA, [10, 23], axis=1)
                matrixA = np.delete(matrixA, [10, 23], axis=0)
                
                # Add time derivative
                difference = np.array([
                    [(g.M / (2*A))], [(g.M / (2*A))], [(b.M / (2*A))], [(b.M / (2*A))],
                    [(c.M / (2*A))], [(c.M / (2*A))], [(d.M / (2*A))], [(d.M / (2*A))],
                    [(e.M / (2*A))], [(e.M / (2*A))], [(f.M / (2*A))], [(f.M / (2*A))],
                ])
                
                for j in range(len(difference)):
                    matrixA[j][j] += difference[j][0] / TIME_STEP
                    matrixA[j][j+12] -= difference[j][0] / TIME_STEP
                    matrixA[j+12][j] += difference[j][0] / TIME_STEP
                    matrixA[j+12][j+12] -= difference[j][0] / TIME_STEP
                
                # Solve
                T = np.dot(np.linalg.inv(matrixA), matrixB)
                
                Teg.append(T[0][0]); Tig.append(T[1][0]); Teb.append(T[2][0]); Tib.append(T[3][0])
                Tec.append(T[4][0]); Tic.append(T[5][0]); Ted.append(T[6][0]); Tid.append(T[7][0])
                Tee.append(T[8][0]); Tie.append(T[9][0]); Tp.append(Tpchange)
                Tef.append(T[10][0]); Tif.append(T[11][0])
                
                Teg.append(T[12][0]); Tig.append(T[13][0]); Teb.append(T[14][0]); Tib.append(T[15][0])
                Tec.append(T[16][0]); Tic.append(T[17][0]); Ted.append(T[18][0]); Tid.append(T[19][0])
                Tee.append(T[20][0]); Tie.append(T[21][0]); Tp.append(Tpchange)
                Tef.append(T[22][0]); Tif.append(T[23][0])
                
                # Calculate total heat
                totalHeat += ((e.kperd * (Tee[2*i] - Tie[2*i]) + h['h_radep'] * (Tie[2*i] - Tp[2*i])
                              + f.kperd * (Tef[2*i] - Tif[2*i]) - h['h_radpf'] * (Tp[2*i] - Tef[2*i]))
                             + (e.kperd * (Tee[2*i+1] - Tie[2*i+1]) + h['h_radep'] * (Tie[2*i+1] - Tp[2*i+1])
                              + f.kperd * (Tef[2*i+1] - Tif[2*i+1]) - h['h_radpf'] * (Tp[2*i+1] - Tef[2*i+1]))) * 60
                
                i += 1
                if i >= 380:
                    break
            
            if totalHeat > 0:
                Tp[2*i-1] = Tpchange - 0.5
                if s == 1:
                    state[0] = 1 - state[0]
                else:
                    state[1] = 1 - state[1]
                print('change')
            else:
                Tp[2*i-1] = Tpchange - 0.5
                print('not change')
    
    # -------------------------------------------------------------------------
    # Post-Processing
    # -------------------------------------------------------------------------
    Tc = [(Tec[i] + Tic[i]) / 2 for i in range(len(Tec))]
    nTc = np.array(Tc)
    Tcs = smooth(nTc - 300, 35) + 300
    
    P = power(Tc, G)
    Ps = power(Tcs, G)
    E = [np.trapz(P[:i]) for i in range(len(P))]
    
    P_increase = [Ps[i] - Ps_NoPCM[i] for i in range(len(P))]
    P_increase_s = smooth(P_increase, 100)
    T_increase = [Tc[i] - Tc_NoPCM[i] for i in range(len(Tc))]
    
    # Print results
    ax = plt.subplot()
    print(np.mean(P[330:750]))
    print(np.mean(P_NoPCM[330:750]))
    print(np.max(P_increase[330:750]))
    print(np.max(P_NoPCM[330:750]))
    print(h['h_radgb'])
    
    # Store optimization results
    P_avg.append(np.mean(Ps[330:750]) - 85.7348512)
    Cost.append((p1.price * r) + (p2.price * (1 - r)))
    Cost_eff.append((np.mean(Ps[330:750]) - 85.7348512) / 85.7348512 * 100 / ((p1.price * r) + (p2.price * (1 - r))))


# =============================================================================
# Final Results
# =============================================================================
print('')
print(np.max(P_avg))
print(P_avg.index(np.max(P_avg)))
print(np.max(Cost_eff))
print(Cost_eff.index(np.max(Cost_eff)))

plt.plot(ratio, Cost_eff)
plt.show()
plt.plot(ratio, P_avg)
plt.show()

print(Tcs[330:750].tolist())
print(Ps[330:750])
