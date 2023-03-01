import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def R(t):
    return t / 100

def G_dt_pend(alpha_G, L, B_l, n_1, Lm1, B_l1, n_2, y_g, GFP):
    dydt = alpha_G / (1 + (L / B_l) ** n_1) + alpha_G / (1 + (Lm1 / B_l1) ** n_2) - y_g * GFP
    return dydt

def L_dt_pend(alpha_L, C, B_c, n_3, y_l, L):
    dydt = alpha_L / (1 + (C / B_c) ** n_3) - y_l * L
    return dydt

def C_dt_pend(alpha_C, Q_r, R, n_4, y_c, C):
    dydt = alpha_C / (1 + (Q_r / R) ** n_4) - y_c * C
    return dydt

def Lm1_dt_pend(alpha_Lm1, Q_r, R, n_5, y_l, Lm1):
    dydt = alpha_Lm1 / (1 + (Q_r / R) ** n_5) - y_l * Lm1
    return dydt

def equation_integration(y, t, alpha_G, B_l, n_1, B_l1, n_2, y_c,
                         alpha_L, B_c, n_3, y_l, alpha_C,
                         Q_r, R, n_4, alpha_Lm1, n_5, y_g):

    GFP, L, C, Lm1 = y
    return [G_dt_pend(alpha_G, L, B_l, n_1, Lm1, B_l1, n_2, y_g, GFP),
            L_dt_pend(alpha_L, C, B_c, n_3, y_l , L),
            C_dt_pend(alpha_C, Q_r, R, n_4, y_c , C),
            Lm1_dt_pend(alpha_Lm1, Q_r, R, n_5, y_l , Lm1)]

t = np.linspace(0, 10, 101)
GFP = 4.231
L = 1.0000000234241
C = 6.443452
Lm1 = 1.022313

y0 = [GFP, L, C, Lm1] 
Q_r = 0.0232312
alpha_G = 1
B_l = 15.03451
n_1 = 5.34501
B_l1 = 3.34501
n_2 = 12.5601
y_c = 0.54301
alpha_L = 0.12301
B_c = 2.23434501
n_3 = 4
alpha_C = 1.01231
Q_r = 0.5
R =0.5
n_4 = 1.04241
y_l =0
n_5 = 1.12311
alpha_Lm1 = 2.0232121
y_g = 1.1231231
sol = odeint(equation_integration, y0, t, args=(alpha_G, B_l, n_1, B_l1, n_2, y_c,
                         alpha_L, B_c, n_3, y_l, alpha_C,
                         Q_r, R, n_4, alpha_Lm1, n_5, y_g))


plt.plot(t, sol[:, 0], 'b', label = 'GFP')
plt.plot(t, sol[:, 1], 'g', label = 'LacI')
plt.plot(t, sol[:, 2], 'r', label = 'CI')
plt.plot(t, sol[:, 3], 'y', label = 'LacIm1')
plt.legend(bbox_to_anchor = (0.6, 1))
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.grid()
plt.show()
