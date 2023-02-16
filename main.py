import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def R(t):
    return t / 100


def G_dt_pend(alpha_G, LacI, B_l, n_1, LacIm1, B_l1, n_2, y_g, GFP):
    dydt = alpha_G / (1 + (LacI / B_l) ** n_1) + alpha_G / (1 + (LacIm1 / B_l1) ** n_2) - y_g * GFP
    return dydt

def L_dt_pend(alpha_L, C_I, B_c, n_3, y_l, L):
    dydt = alpha_L / (1 + (C_I / B_c) ** n_3) - y_l * L
    return dydt

def C_dt_pend(alpha_C, Q_r, R, n_4, y_c, C):
    dydt = alpha_C / (1 + (Q_r + R) **  n_4) - y_c * C
    return dydt

def Lm1_dt_pend(alpha_Lm1, Q_r, R, n_5, y_lm1, M1):
    dydt = alpha_Lm1 / (1 + (Q_r + R) ** n_5) - y_lm1 * m1
    return dydt


def equation_integration(y, t, alpha_G, B_l, n_1, B_l1, n_2, y_c,
                         alpha_L, B_c, n_3, y_l, L, alpha_C,
                         Q_r, R, n_4, C, alpha_Lm1, n_5, m1):
    GFP, LacI, C_I, LacIm1 = y
    return [G_dt_pend(alpha_G, LacI, B_l, n_1, LacIm1, B_l1, n_2, y_c, GFP),
            L_dt_pend(alpha_L, C_I, B_c, n_3, y_l , L),
            C_dt_pend(alpha_C, Q_r, R, n_4, y_c , C),
            Lm1_dt_pend(alpha_Lm1, Q_r, R, n_5, y_l , M1)]

t = np.linspace(0, 10, 101)
y0 = [np.pi - 0.23421, 0.234241, 0.3453452, 0.022313]

Q_r = 0.0232312
alpha_G = 0.34501
B_l = 0.03451
n_1 = 0.34501
B_l1 = 0.34501
n_2 = 0.5601
y_c = 0.54301
alpha_L = 0.12301
B_c = 0.23434501
n_3 = 0.3543401
alpha_C = 0.01231
Q_r = 0.12301
R = 0.0123131
n_4 = 0.04241
y_l = 0.0422
n_5 = 0.012311
L = 0.33335435
C = 0.09121
alpha_Lm1 = 0.232121
m1 = 0.1212313

sol = odeint(equation_integration, y0, t, args=(alpha_G, B_l, n_1, B_l1, n_2, y_c,
                         alpha_L, B_c, n_3, y_l, L, alpha_C,
                         Q_r, R, n_4, C, alpha_Lm1, n_5, m1))



plt.plot(t, sol[:, 0], 'b', label = 'LacIm1 and LacI repress GFP')
plt.plot(t, sol[:, 1], 'g', label = 'CI represses LacI')
plt.plot(t, sol[:, 2], 'r', label = 'LuxR activates the CI gene')
plt.plot(t, sol[:, 3], 'y', label = 'LuxR activates the LacIm1 gene')
plt.legend(loc = 'best')
plt.xlabel('t')
plt.ylabel('f(x)')
plt.grid()
plt.show()
