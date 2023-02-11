from scipy.integrate import odeint
import np

def R(t):
    return t/100


def G_dt_pend(y):
    alpha_g, LacI, B_l,n_1,LacIm1,B_l1,n_2,y_c, GFP = y
    dydt = alpha_g/(1+(LacI/B_l)**n_1) + alpha_g/(1+(LacIm1/B_l1)**n_2) - y_g*GFP
    return dydt

def L_dt_pend(y):
    alpha_l, C_I, B_c,n_3, y_l, L = y
    dydt = alpha_l/(1+(C_I/B_c)**n_3) - y_l*L
    return dydt

def C_dt_pend(y):
    alpha_c, Q_r, R,n_4, y_c, C = y
    dydt = alpha_c/(1+(Q_r+R)**n_4) - y_c*C
    return dydt

def Lm1_dt_pend(y):
    alpha_lM1, Q_r, R,n_5,y_lm1, M1 = y
    dydt = alpha_lM1/(1+(Q_r+R)**n_5) - y_lm1*M1
    return dydt


def equation_integration(y, t, alpha_g, B_l,n_1,B_l1,n_2,y_c,
                         alpha_l, B_c,n_3,y_lL, alpha_c,
                         Q_r, R,n_4, y_cC, alpha_lM1,n_5,y_lm1M1):
    GFP, LacI, C_I, LacIm1 = y
    return [G_dt_pend([alpha_g, LacI, B_l,n_1,LacIm1,B_l1,n_2,y_c, GFP]),
            L_dt_pend([alpha_l, C_I, B_c,n_3, y_lL]),
            C_dt_pend([alpha_c, Q_r, R,n_4, y_cC]),
            Lm1_dt_pend([alpha_lM1, Q_r, R,n_5,y_lm1M1])]

t = np.linspace(0, 10, 101) 
y0 = [np.pi - 0.23421, 0.234241,0.3453452,0.022313] 
Q_r = 0.0232312
alpha_g= 0.34501
B_l= 0.03451
n_1= 0.34501
B_l1= 0.34501
n_2= 0.5601
y_c= 0.54301
alpha_l= 0.12301
B_c= 0.23434501
n_3= 0.3543401
y_lL= 0.12301
alpha_c= 0.01231
Q_r= 0.12301
R= 0.0123131
n_4= 0.04241
y_cC= 0.0234241
alpha_lM1= 0.0234211
n_5= 0.012311
y_lm1M1= 0.05551

t = np.linspace(0, 10, 101)
y0 = [np.pi - 0.1, 0.0,0,0]
sol = odeint(equation_integration, y0, t, args=(alpha_g, 
                         B_l,n_1,B_l1,n_2,y_c,
                         alpha_l, B_c,n_3,y_lL, alpha_c,
                         Q_r, R,n_4, y_cC, alpha_lM1,n_5,y_lm1M1))


import matplotlib.pyplot as plt
plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
