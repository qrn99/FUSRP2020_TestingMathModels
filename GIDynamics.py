"""
Class of functions simulating Glucose and Insuline dynamics
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

k_js = 10**(-0.8)
k_gj = 10**(-1)
k_jl = 10**(-1.5)
tao = 10**2
k_gl = 10**(-0.8)
k_xg = 10**(-1.9)
k_xgi = 10**(-7)
eta = 10**(-1.9)
k_lambda = 10**(-1.3)
f_gj = 10**0.8
k_xi = 10**(-1.9)
gamma = 10**0.2
beta = 10**2
D = 10**2.5
G_b = 10
I_b = 10**1.8
G_prod0 = -G_b

tau = 50

def gi(x, t):
    """
    gi dynamics complete
    :param x:
    :type x:
    :param t:
    :type t:
    :return:
    :rtype:
    """
    s = x[0]
    j = x[1]
    l = x[2]
    g = x[3]
    i = x[4]
    gprod = k_lambda/(k_lambda/G_prod0+(g-G_b))

    phi = 0 if t < tao else j * (t - tau)

    dsdt = -k_js * s
    djdt = k_js * s - k_gj * j - k_jl * j
    dldt = k_jl * phi - k_gl * l
    dgdt = -(k_xg + k_xgi*i)*g + gprod + eta*(k_gj*j + k_gl*l)
    didt = k_xi*I_b*((beta**gamma+1)/(beta**gamma*(G_b/(g+f_gj*(k_gj*j + k_gl*l)))**gamma + 1)
                     - i/I_b)
    return [dsdt, djdt, dldt, dgdt, didt]

x0 = [70, 0, 0, G_b, I_b]
t = np.linspace(0, 2000, 5000)
x = odeint(gi, x0, t)

s = x[:, 0]
j = x[:, 1]
l = x[:, 2]
g = x[:, 3]
i = x[:, 4]

plt.plot(t, s, label='S')
plt.plot(t, j, label='J')
plt.plot(t, l, label='L')
plt.plot(t, g, label='G')
plt.plot(t, i, label='I')
plt.legend()
plt.show()
