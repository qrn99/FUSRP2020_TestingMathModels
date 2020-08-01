"""
Class of functions simulating Glucose and Insuline dynamics
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from typing import List

k_js = 0.1
k_gj = 0.1
k_jl = 0.0316
tau = 50
k_gl = 0.1
k_xg = 0.01
k_xgi = 10**(-6.5)
eta = 10**(-1.9)
k_lambda = 0.0316
f_gj = 10
k_xi = 0.01
gamma = 5
beta = 100
D = 10**2.5
G_b = 7
I_b = 1.7*10**(-7)
G_prod0 = k_xgi*G_b*I_b+k_xg*G_b
J = [0]*120


def time_delay(t, J):
    if t < tau:
        return 0
    else:
        return J[round(t-tau)]


def Gprod(g):
    return k_lambda / (k_lambda / G_prod0 + (g - G_b))


def gi(x, t):
    """
    gi dynamics complete
    args are the input params
    """
    # k_js, k_gj, k_jl, tau, k_gl, k_xg, k_xgi, eta, k_lambda, f_gj, k_xi, gamma, beta, D, G_b, I_b = args

    s = x[0]
    j = x[1]
    l = x[2]
    g = x[3]
    i = x[4]

    dsdt = -k_js * s
    djdt = k_js * s - k_gj * j - k_jl * j
    J[round(t)] = j
    dldt = k_jl * time_delay(t, J) - k_gl * l
    dgdt = -(k_xg + k_xgi * i)*g + Gprod(g) + eta*(k_gj*j + k_gl*l)
    num = beta**gamma + 1
    G_tilda = g + f_gj * (k_gj * j + k_gl * l)
    denom = (beta**gamma) * (G_b/G_tilda)**gamma + 1
    didt = k_xi * I_b * (num/denom - i/I_b)
    return [dsdt, djdt, dldt, dgdt, didt]


x0 = [417, 0, 0, G_b, I_b]
t = np.linspace(0, 119, 119)
x = odeint(gi, x0, t)

s = x[:, 0]
j = x[:, 1]
l = x[:, 2]
g = x[:, 3]
i = x[:, 4]

# plt.plot(t, s, label='S')
# plt.plot(t, j, label='J')
# plt.plot(t, l, label='L')
plt.plot(t, g, label='G')
# plt.plot(t, i, label='I')
plt.legend()
plt.show()
