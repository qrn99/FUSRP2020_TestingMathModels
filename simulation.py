"""
The main simulation class for fitted glucose curve (and other curves)
"""
import pandas as pd
from typing import List, Tuple
import GIDynamics
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from itertools import product
import Fitting_Params

if __name__ == '__main__':
    g_exp = pd.read_csv('file_name')

    k_js = None
    k_gj = None
    k_jl = None
    tau = None
    k_gl = None
    k_xg = None
    k_xgi = None
    eta = None
    k_lambda = None
    f_gj = None
    k_xi = None
    gamma = None
    beta = None
    D = None
    G_b = None
    I_b = None
    c = None

    ub_lb = {k_js: [10**(-0.9), 10**(-0.7)],
             k_gj: [10**(-1.5), 10**(-0.5)],
             k_jl: [10**(-2), 10**(-1)],
             tau: [50, 150],
             k_gl: [10**(-0.9), 10**(-0.7)],
             k_xg: [10**(-2.0), 10**(-1.8)],
             k_xgi: [10**(-8), 10**(-7)],
             eta: [10**(-1.8), 10**(-2)],
             k_lambda: [10**(-1.4), 10**(-1.2)],
             f_gj: [10**0.7, 10**0.9],
             k_xi: [10**(-1.8), 10**(-2)],
             gamma: [10**0.1, 10**0.3],
             beta: [50, 100],
             D: [100, 500],
             G_b: [10**0.7, 10],
             I_b: [10**1.7, 10**1.9],
             # G_prod0: [None, None]
             }

    theta0 = Fitting_Params.get_theta0(ub_lb)

    theta_pool = Fitting_Params.get_theta_pool(theta0, 5)

    params, loss = Fitting_Params.fit_params(theta0, theta_pool, g_exp)

    # Produce the glucose simulation below
    x0 = [D, 0, 0, G_b, I_b] # initial condition
    t = np.linspace(0, len(g_exp), len(g_exp))
    x = odeint(GIDynamics.g, x0, t)

    s = x[:, 0]
    j = x[:, 1]
    l = x[:, 2]
    g = x[:, 3]
    i = x[:, 4]

    plt.plot(t, g, label='G')
    plt.legend()
    plt.show()

    # Show all the curves
    plt.plot(t, s, label='S')
    plt.plot(t, j, label='J')
    plt.plot(t, l, label='L')
    plt.plot(t, g, label='G')
    plt.plot(t, i, label='I')

plt.legend()
plt.show()
