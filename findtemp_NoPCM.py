
# * Concept of this method
# * * dT/dt = (T[i+1]-T[i])/120
# * * When A*X = B, X = invA*B

# * -- Import libraries and data

import numpy as np
import matplotlib.pyplot as plt
from data_smooth import Tas, Gs
from power import power

# * -- Class of properties


class Material:

    def properties(self, c_p, d, emi, p, k):
        A = 1.61
        self.c_p = c_p
        self.d = d
        self.emi = emi
        self.m = p*A*d
        self.k = k
        self.M = self.m*c_p
        self.kperd = k/d

# * -- Material's propeties


g = Material()
b = Material()
c = Material()
d = Material()
e = Material()
f = Material()

g.properties(c_p=840, d=0.0032, emi=0.85, p=2520, k=0.96)
b.properties(c_p=3135, d=0.0005, emi=0.53, p=934, k=0.24)
c.properties(c_p=710.08, d=0.0003, emi=0.75, p=2330, k=148)
d.properties(c_p=3135, d=0.0005, emi=0.53, p=934, k=0.24)
e.properties(c_p=1050, d=0.0005, emi=0.08, p=2700, k=0.1583)
f.properties(c_p=840, d=0.005, emi=0.85, p=2500, k=0.96)

A = 1.61
h_cong = 24.486272
h_conf = 24.486272
state = 0

# * Important parameters
sigma = 5.67 * 10**(-8)
G = Gs/116
Pg = 0.925 * 0.9 * G
T_sky = []
for item in Tas:
    T_sky.append(0.0552*(item**1.5))

chagestate = False

# * -- Lists of temperature through time

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
Tef = [Tas[0], Tas[1]]
Tif = [Tas[0], Tas[1]]

# * -- Calculate in each 2 timepoints
# * 1500 timepoint 2 timepoints per time = 786 times

i = 1
while i < 751:

    # * -- H radiation

    h_radag = sigma * g.emi * \
        (Teg[2*i-1]**2 + T_sky[2*i-1]**2) * (Teg[2*i-1] + T_sky[2*i-1])
    h_radgb = (sigma * (Teb[2*i-1]**2 + Tig[2*i-1]**2)
               * (Teb[2*i-1] + Tig[2*i-1]))/((1/g.emi)+(1/b.emi)-1)
    h_radbc = (sigma * (Tec[2*i-1]**2 + Tib[2*i-1]**2)
               * (Tec[2*i-1] + Tib[2*i-1]))/((1/b.emi)+(1/c.emi)-1)
    h_radcd = (sigma * (Ted[2*i-1]**2 + Tic[2*i-1]**2)
               * (Ted[2*i-1] + Tic[2*i-1]))/((1/c.emi)+(1/d.emi)-1)
    h_radde = (sigma * (Tee[2*i-1]**2 + Tid[2*i-1]**2)
               * (Tee[2*i-1] + Tid[2*i-1]))/((1/d.emi)+(1/e.emi)-1)
    h_radef = (sigma * (Tef[2*i-1]**2 + Tie[2*i-1]**2)
               * (Tef[2*i-1] + Tie[2*i-1]))/((1/e.emi)+(1/f.emi)-1)
    h_radfa = sigma * f.emi * \
        (Tas[2*i-1]**2 + Tif[2*i-1]**2) * (Tas[2*i-1] + Tif[2*i-1])

    h_radg = (sigma * (Teg[2*i-1]**2 + Tig[2*i-1]**2) *
              (Teg[2*i-1] + Tig[2*i-1]))/((1/g.emi)+(1/g.emi)-1)
    h_radb = (sigma * (Teb[2*i-1]**2 + Tib[2*i-1]**2) *
              (Teb[2*i-1] + Tib[2*i-1]))/((1/b.emi)+(1/b.emi)-1)
    h_radc = (sigma * (Tec[2*i-1]**2 + Tic[2*i-1]**2) *
              (Tec[2*i-1] + Tic[2*i-1]))/((1/c.emi)+(1/c.emi)-1)
    h_radd = (sigma * (Ted[2*i-1]**2 + Tid[2*i-1]**2) *
              (Ted[2*i-1] + Tid[2*i-1]))/((1/d.emi)+(1/d.emi)-1)
    h_rade = (sigma * (Tee[2*i-1]**2 + Tie[2*i-1]**2) *
              (Tee[2*i-1] + Tie[2*i-1]))/((1/e.emi)+(1/e.emi)-1)
    h_radf = (sigma * (Tef[2*i-1]**2 + Tif[2*i-1]**2) *
              (Tef[2*i-1] + Tif[2*i-1]))/((1/f.emi)+(1/f.emi)-1)

    # * -- matrix B construction

    B0_1 = [-(Pg[2*i]/2)-(h_radag*T_sky[2*i])-(h_cong*Tas[2*i])]
    B0_2 = [-(Pg[2*i+1]/2)-(h_radag*T_sky[2*i+1])-(h_cong*Tas[2*i+1])]
    B1_1 = [-(Pg[2*i]/2)]
    B1_2 = [-(Pg[2*i+1]/2)]
    B2, B3, B4, B5, B6, B7, B8, B9, B10 = [0], [
        0], [0], [0], [0], [0], [0], [0], [0]
    B11_1 = [-(h_radfa*Tas[2*i])-(h_conf*Tas[2*i])]
    B11_2 = [-(h_radfa*Tas[2*i+1])-(h_conf*Tas[2*i+1])]

    matrixB = np.array([B0_1, B1_1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11_1,
                       B0_2, B1_2, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11_2])
    
    # * k between layer
    kperdgb = 1/((1/g.kperd)+(1/b.kperd))
    kperdbc = 1/((1/b.kperd)+(1/c.kperd))
    kperdcd = 1/((1/c.kperd)+(1/d.kperd))
    kperdde = 1/((1/d.kperd)+(1/e.kperd))
    kperdef = 1/((1/e.kperd)+(1/f.kperd))

    # * -- matrix A construction
    # * Equation will be solved 2 time points in each times

    A0 = [-g.kperd-h_radag-h_cong-h_radg, g.kperd +
            h_radg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    A1 = [-kperdgb+2*g.kperd+h_radg, -2*g.kperd-h_radg-h_radgb, h_radgb +
            b.kperd, kperdgb-b.kperd, 0, 0, 0, 0, 0, 0, 0, 0]
    A2 = [kperdgb-g.kperd, g.kperd+h_radgb, -b.kperd-h_radgb-b.kperd -
            h_radb, -kperdgb+b.kperd+b.kperd+h_radb, 0, 0, 0, 0, 0, 0, 0, 0]
    A3 = [0, 0,-kperdbc+2*b.kperd+h_radb, -2*b.kperd-h_radb-h_radbc, h_radbc +
            c.kperd, kperdbc-c.kperd, 0, 0, 0, 0, 0, 0]
    A4 = [0, 0, kperdbc-b.kperd, b.kperd+h_radbc, -c.kperd-h_radbc -
            c.kperd - h_radc, -kperdbc+c.kperd+c.kperd+h_radc, 0, 0, 0, 0, 0, 0]
    A5 = [0, 0, 0, 0,-kperdcd+2*c.kperd+ h_radc, -2*c.kperd-h_radc-h_radcd,
            h_radcd + d.kperd, kperdcd-d.kperd, 0, 0, 0, 0]
    A6 = [0, 0, 0, 0, kperdcd-c.kperd, c.kperd+h_radcd, -d.kperd-h_radcd -
            d.kperd - h_radd, -kperdcd+d.kperd+d.kperd+h_radd, 0, 0, 0, 0]
    A7 = [0, 0, 0, 0, 0, 0,-kperdde+2*d.kperd+ h_radd, -2*d.kperd-h_radd -
            h_radde, h_radde + e.kperd, kperdde-e.kperd, 0, 0]
    A8 = [0, 0, 0, 0, 0, 0, kperdde-d.kperd, d.kperd+h_radde, -e.kperd -
            h_radde-e.kperd - h_rade, -kperdde+e.kperd+e.kperd+h_rade, 0, 0]
    A9 = [0, 0, 0, 0, 0, 0, 0, 0,-kperdef+2*e.kperd+ h_rade, -2*e.kperd - 
          h_rade-h_radef, h_radef + f.kperd, kperdef-f.kperd]
    A10 = [0, 0, 0, 0, 0, 0, 0, 0, kperdef-e.kperd, e.kperd+h_radef, -
           f.kperd - h_radef-f.kperd-h_radf, -kperdef+f.kperd+f.kperd+h_radf]
    A11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, f.kperd +
           h_radf-h_conf, -f.kperd-h_radf-h_radfa]

    vectorOf0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    matrixA = np.array([np.append(A0, vectorOf0),
                        np.append(A1, vectorOf0),
                        np.append(A2, vectorOf0),
                        np.append(A3, vectorOf0),
                        np.append(A4, vectorOf0),
                        np.append(A5, vectorOf0),
                        np.append(A6, vectorOf0),
                        np.append(A7, vectorOf0),
                        np.append(A8, vectorOf0),
                        np.append(A9, vectorOf0),
                        np.append(A10, vectorOf0),
                        np.append(A11, vectorOf0),

                        np.append(vectorOf0, A0),
                        np.append(vectorOf0, A1),
                        np.append(vectorOf0, A2),
                        np.append(vectorOf0, A3),
                        np.append(vectorOf0, A4),
                        np.append(vectorOf0, A5),
                        np.append(vectorOf0, A6),
                        np.append(vectorOf0, A7),
                        np.append(vectorOf0, A8),
                        np.append(vectorOf0, A9),
                        np.append(vectorOf0, A10),
                        np.append(vectorOf0, A11),
                        ])

    # * -- add dT/dt or (T[i+1]-T[i])/120 to matrix A

    diference = np.array([
        [((g.M)/(2*A))],
        [((g.M)/(2*A))],
        [(b.M)/(2*A)],
        [(b.M)/(2*A)],
        [(c.M)/(2*A)],
        [(c.M)/(2*A)],
        [(d.M)/(2*A)],
        [(d.M)/(2*A)],
        [(e.M)/(2*A)],
        [(e.M)/(2*A)],
        [(f.M)/(2*A)],
        [((f.M)/(2*A))],
    ])

    for j in range(len(diference)):
        matrixA[j][j] = matrixA[j][j]-diference[j][0]/120
        matrixA[j][j+12] = matrixA[j][j+12]+diference[j][0]/120
        matrixA[j+12][j] = matrixA[j+12][j]-diference[j][0]/120
        matrixA[j+12][j+12] = matrixA[j+12][j+12]+diference[j][0]/120

    # * -- Get vector X

    matrixinvA = np.linalg.inv(matrixA)
    T = np.dot(matrixinvA, matrixB)
    Teg.append(T[0][0])
    Tig.append(T[1][0])
    Teb.append(T[2][0])
    Tib.append(T[3][0])
    Tec.append(T[4][0])
    Tic.append(T[5][0])
    Ted.append(T[6][0])
    Tid.append(T[7][0])
    Tee.append(T[8][0])
    Tie.append(T[9][0])
    Tef.append(T[10][0])
    Tif.append(T[11][0])
    Teg.append(T[12][0])
    Tig.append(T[13][0])
    Teb.append(T[14][0])
    Tib.append(T[15][0])
    Tec.append(T[16][0])
    Tic.append(T[17][0])
    Ted.append(T[18][0])
    Tid.append(T[19][0])
    Tee.append(T[20][0])
    Tie.append(T[21][0])
    Tef.append(T[22][0])
    Tif.append(T[23][0])

    i = i+1


# * -- Find Temperature of solar panel

Tc_NoPCM = []
for i in range(len(Tec)):
    Tc_NoPCM.append((Tec[i]+Tic[i])/2)


def smooth(x, box_pts):
    box = np.ones(box_pts)/box_pts
    x_smooth = np.convolve(x, box, mode='same')
    return x_smooth


nTc = np.array(Tc_NoPCM)
Tcs_NoPCM = smooth(nTc - 300, 35) + 300

# * -- Find power and total energy of solar panel
P_NoPCM = power(Tc_NoPCM, G)
Ps_NoPCM = power(Tcs_NoPCM, G)
E_NoPCM = []
for i in range(len(P_NoPCM)):
    E_NoPCM.append(np.trapz(P_NoPCM[:i]))

# * -- Data analyse
print(Gs[330:750].tolist())
# x = range(1500)
# ax = plt.subplot()
# ax2 = plt.twinx(ax)
# ax2.plot(x, Tas[:-2])
# ax2.plot(x, Teg[:-2])
# ax2.plot(x, Teb[:-2])
# ax2.plot(x, Tec[:-2])
# ax2.plot(x, Ted[:-2])
# ax2.plot(x, Tee[:-2])
# ax2.plot(x, Tef[:-2])
# ax2.legend(['a','g','b','c','d','e','f'])
# print(h_radgb)
# plt.show()
