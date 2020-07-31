"""
Class of functions that fits the best parameters to the system of diff eqns
"""
import pandas as pd
from typing import List, Tuple
import GIDynamics
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from itertools import product

# import the whole glucose data file
data = pd.read_csv(file='data_clean.csv')
data.head()


# Determine the time_scope to examine, default is 120min (in the paper)
time_scope = 120
# find 120 min of glucose data to study
# set the first time point for the slice of data to be 0
# rename each time point as consecutive integers starting from 0

m = [None]*time_scope # time points going from 0 to 120
print(len(m)) # check whether m is of length 120, if not, there there are data missing for some minutes

G_data = [] # enter the data per minute, if there is data missing, use interpolation


def interpolation(lb, ub, num_pt):
    # TODO: implement the function
    return []


def unpack(x: list):
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = x
    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16


# the complete glucose data in 120 mins to use
G = [100]*120

# Initialize arguments
k_js = None
k_gj = None
k_jl = None
tao = None
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
G_prod0 = None

units_dict = {}

# set up upper and lower bounds for the parameters in a dict
ub_lb = {k_js: [None, None],
         k_gj: [None, None],
         k_jl: [None, None],
         tao: [None, None],
         k_gl: [None, None],
         k_xg: [None, None],
         k_xgi: [None, None],
         eta: [None, None],
         k_lambda: [None, None],
         f_gj: [None, None],
         k_xi: [None, None],
         gamma: [None, None],
         beta: [None, None],
         D: [None, None],
         G_b: [None, None],
         I_b: [None, None],
         # G_prod0: [None, None]
         }


def compute_loss(x: List[float], g: List[float]) -> float:
    """
    Formula:
    J_exp(theta, alpha) = 1/(5*(G_max-G_min)^2)*(sum((G_exp(i)-G_num(Ti))^2)) + error of I
    """
    m = len(x)
    g_max = max(g)
    g_min = min(g)
    loss = 0

    for i in range(m):
        loss += (g[i] - x[i]) ** 2

    return 1/(m*(g_max-g_min)**2)*loss


def model_loss(params: List[float], g: List[float]) -> float:
    """
    g is the experimental glucose values
    """
    x0 = params
    num_points = len(g)
    t = np.linspace(0, num_points, num_points)
    x = odeint(GIDynamics.gi, x0, t)
    fitted_glucose = list(x[:, 3])
    loss = compute_loss(fitted_glucose, g)
    return loss


def get_theta0(bounds_dict: dict) -> List[float]:
    """
    return the initial guess for theta by taking the middle point of the possible range
    """
    return [sum(bounds_dict[x])/len(ub_lb[x]) for x in bounds_dict.keys()]


def get_theta_pool(theta: List[float], size) -> List[List[float]]:
    """
    Return the a pool of values for each params in theta to sample from.
    There are size number of candidate for each param.
    """
    # TODO: can consider using unit_dict to keep track of how much to adjust for each param.
    return [interpolation(x-1, x+1, size) for x in theta]


def fit_params(theta0, theta_pool, g) -> List[List[float], float]:
    """
    return the list of optimal params and the minimized loss value
    """
    min_loss = model_loss(theta0, g)
    best_param = None

    for p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 \
            in product(unpack(theta_pool)):
        params = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16]
        loss = model_loss(params, G)
        if loss < min_loss:
            min_loss = loss
            best_param = params
    return [best_param, min_loss]

# TODO: do something to counteract measurement error

#
# def hypothesis(init_cond: List[float], length: int=120) -> List:
#     x0 = init_cond
#     t = np.linspace(0, length, length)
#     x = odeint(GIDynamics.gi, x0, t)
#     return list(x[:, 3])
#
#
# def gradient_descent(theta: List[float], num_iters: int, lr: float=1e5):
#     """
#     Implement gradient descent to minimize the losss
#     :return:
#     :rtype:
#     """
#     cost = np.ones(num_iters)
#     for i in range(0, num_iters):
#         theta[0]
#
# def get_grad_vec(g_exp) -> List:
#     """
#     Return a list of partial derivatives wrt the params.
#     :param g_exp:
#     :type g_exp:
#     :return:
#     :rtype:
