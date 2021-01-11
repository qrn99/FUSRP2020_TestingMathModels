"""
Class of functions that fits the best parameters to the system of diff eqns
"""
from datetime import datetime
import pandas as pd
from typing import List, Tuple, Dict, Any, Union
import GIDynamics
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from itertools import product
import math

# import the whole glucose data file
data = pd.read_csv('2020-07-2021_lunch_peak.csv')
data.head()

print(type(data.time))
print(type(data.glucose))

# date_time_obj = datetime. strptime(date_time_str, '%d/%m/%y %H:%M:%S')

# Convert the time str to datetime object
time_list = [datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in data.time]

# Replace the time str using the new datetime object in dataframe
data.time = time_list

# Determine the time_scope to examine, default is 120min (in the paper)
time_scope = data.shape[0]  # 173 right now

g_exp = list(data.glucose)

print(len(g_exp) == 173)

# Initialize arguments. This is the keys they will use to find their initialized values in initial dict
k_js = 0
k_gj = 1
k_jl = 2
tau = 3
k_gl = 4
k_xg = 5
k_xgi = 6
eta = 7
k_lambda = 8
f_gj = 9
k_xi = 10
gamma = 11
beta = 12
D = 13
G_b = 14
I_b = 15
G_prod0 = 16
J = [0] * len(g_exp)

# set up initial values for the parameters: remember to replace the list with values
init_dict = {0: 0.1,
             1: 0.1,
             2: 0.0316,
             3: 100,
             4: 0.1,
             5: 0.01,
             6: 10 ** (-6.5),
             7: 10 ** (-1.9),
             8: 0.0316,
             9: 10,
             10: 0.01,
             11: 5,
             12: 100,
             13: 10 ** 3,
             14: 10 ** 0.8,
             15: 10 ** (1.7 - 9)
             # G_prod0: [None, None]
             }


def unpack(x: Dict):
    y = list(x.values())
    p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15 = y
    return p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15

# the complete glucose data in 120 mins to use
# G = [100]*120

# Fit curves and output fitted glucose
def fit_gi(x, t):
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

    dsdt = -params[0] * s
    djdt = params[0] * s - params[1] * j - params[2] * j
    J[round(t)] = j
    dldt = params[2] * time_delay(params, t, J) - params[4] * l
    dgdt = -(params[5] + params[6] * i) * g + Gprod(params, g) + params[7] * (params[1] * j + params[4] * l)
    num = params[12] ** params[11] + 1
    G_tilda = g + params[9] * (params[1] * j + params[4] * l)
    denom = (params[12] ** params[11]) * (params[14] / G_tilda) ** params[11] + 1
    didt = params[10] * params[15] * (num / denom - i / params[15])
    return [dsdt, djdt, dldt, dgdt, didt]


def compute_loss(x: List[float], g: List[float]) -> float:
    """
    Implement the formula into a function.
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


def interp_range(l: List[Union[float, int]], num_pt: int) -> List[Union[int, float]]:
    """
    Interpolate the list l evently into a length of num_pt.

    >>> l = [1, 4, 5]
    >>> interp_range(l, 5)
    [1, 2, 3, 4, 5]
    """
    if num_pt > 1:
        max_val = max(l)
        min_val = min(l)
        diff = max_val - min_val
        step = diff/(num_pt-1)

        to_return = []
        cur = min_val
        while cur <= max_val:
            to_return.append(cur)
            cur += step
        return to_return
    else:
        to_return = l
        return to_return


def generate_param_range(init_dict: Dict[int, float], num_pt: int) -> Dict[int, List[float]]:
    """
    Generate a dict of list which indicates the suggested parameter pools to choose parameters from.
    The range should vary depending on the parameter's scale. The scale should have an extra significant digit that spans
    the unsearched area.
    For example, if we start with 10, then we search 0~20. [0, 4, 8, 12, 16, 20]
    If a better result is achieved at 16, then we search 12~20. 8/5 = [10, 12, 14, 16, 18, 20]
    """
    if num_pt > 1:
        scale = [math.log10(x) for x in list(init_dict.values())]
        range = [[x-0.2, x+0.2] for x in scale]
        pool = [interp_range(x, num_pt) for x in range]
    else:
        scale = [math.log10(x) for x in list(init_dict.values())]
        range = [[x] for x in scale]
        pool = [interp_range(x, num_pt) for x in range]
    range_pool = []
    for param_pool in pool:
        range_pool.append([10**x for x in param_pool])

    range_dict = {i:v for i, v in enumerate(range_pool)}
    return range_dict


def model_loss(params: List[float], g: List[float]) -> float:
    """
    Return the loss for a model with fitted parameters
    g is the experimental glucose values
    """
    x0 = [params[13], 0, 0, params[14], params[15]]
    J = [0]*len(g)
    num_points = len(g)
    # print(num_points)


    def model_get_Gprod0(params: List[float]) -> float:
        """
        Return the value of Gprod0 given a list of input params
        """
        return params[6] * params[14] * params[15] + params[5] * params[14]

    def model_Gprod(params: List[float], g: float) -> float:
        """
        Calculate the value of G_prod at some point glucose level = g.
        """
        return params[8] / (params[8] / model_get_Gprod0(params) + (g - params[14]))

    def model_time_delay(params: List[float], t: float, J: List[float]) -> float:
        if t < params[3]:
            return 0
        else:
            return J[round(t - params[3])]


    def model_gi(x, t):
        """
        model gi dynamics complete
        """
        # k_js, k_gj, k_jl, tau, k_gl, k_xg, k_xgi, eta, k_lambda, f_gj, k_xi, gamma, beta, D, G_b, I_b = args

        s = x[0]
        j = x[1]
        l = x[2]
        g = x[3]
        i = x[4]

        dsdt = -params[0] * s
        djdt = params[0] * s - params[1] * j - params[2] * j
        J[int(t)] = j
        dldt = params[2] * model_time_delay(params, t, J) - params[4] * l
        dgdt = -(params[5] + params[6] * i) * g + model_Gprod(params, g) + params[7] * (params[1] * j + params[4] * l)
        num = params[12] ** params[11] + 1
        G_tilda = g + params[9] * (params[1] * j + params[4] * l)
        denom = (params[12] ** params[11]) * (params[14] / G_tilda) ** params[11] + 1
        didt = params[10] * params[15] * (num / denom - i / params[15])
        return [dsdt, djdt, dldt, dgdt, didt]

    time = np.linspace(0, num_points-2, num_points-1)
    func = odeint(model_gi, x0, time)
    fitted_glucose = list(func[:, 3])
    loss = compute_loss(fitted_glucose, g)
    return loss


def get_Gprod0(params: List[float]) -> float:
    """
    Return the value of Gprod0 given a list of input params
    """
    return params[6]*params[14]*params[15]+params[5]*params[14]


def fit_params(init_guess: Dict[int, float], g: List[float], num_pt: int) -> [Dict[int, float], float]:
    """
    return the dict of optimal params together with the minimized loss value
    """
    # Generate the pool from which we sample params
    param_pool = generate_param_range(init_guess, num_pt)

    # append Gprod0 to the initial params
    init_params = list(init_guess.values())
    init_params.append(init_guess[6]*init_guess[14]*init_guess[15]+init_guess[5]*init_guess[14])

    # Initialize the min loss of model by the initial params
    min_loss = model_loss(init_params, g)
    best_param = init_params

    cur_it = 0
    total_it = num_pt**16 if num_pt > 1 else 15

    for p0 in param_pool[0]:
        print(len(param_pool[0]))
        for p1 in param_pool[1]:
            for p2 in param_pool[2]:
                for p3 in param_pool[3]:
                    for p4 in param_pool[4]:
                        for p5 in param_pool[5]:
                            for p6 in param_pool[6]:
                                for p7 in param_pool[7]:
                                    for p8 in param_pool[8]:
                                        for p9 in param_pool[9]:
                                            for p10 in param_pool[10]:
                                                for p11 in param_pool[11]:
                                                    for p12 in param_pool[12]:
                                                        for p13 in param_pool[13]:
                                                            for p14 in param_pool[14]:
                                                                for p15 in param_pool[15]:

                                                                # for p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15 \
                                                                #         in unpacked_params:
                                                                    cur_it += 1
                                                                    params = [p0, p1, p2, p3, p4, p5, p6, p7, p8,
                                                                              p9, p10, p11, p12, p13, p14, p15]
                                                                    params.append(get_Gprod0(params))
                                                                    loss = model_loss(params, g)
                                                                    if loss < min_loss:
                                                                        min_loss = loss
                                                                        best_param = params
                                                                        print(min_loss, "Progress:{0}".format(cur_it/total_it))
                                                                    if cur_it % 1000 == 0:
                                                                        to_save = pd.DataFrame(best_param)
                                                                        f_name = '5_after {0} iterations'.format(cur_it) + '.csv'
                                                                        to_save.to_csv(f_name)
        #                                                             if cur_it / total_it >= 0.01:
        #                                                                 break
        #                                                         if cur_it / total_it >= 0.01:
        #                                                             break
        #                                                     if cur_it / total_it >= 0.01:
        #                                                         break
        #                                                 if cur_it / total_it >= 0.01:
        #                                                     break
        #                                             if cur_it / total_it >= 0.01:
        #                                                 break
        #                                         if cur_it / total_it >= 0.01:
        #                                             break
        #                                     if cur_it / total_it >= 0.01:
        #                                         break
        #                                 if cur_it / total_it >= 0.01:
        #                                     break
        #                             if cur_it / total_it >= 0.01:
        #                                 break
        #                         if cur_it / total_it >= 0.01:
        #                             break
        #                     if cur_it / total_it >= 0.01:
        #                         break
        #                 if cur_it / total_it >= 0.01:
        #                     break
        #             if cur_it / total_it >= 0.01:
        #                 break
        #         if cur_it / total_it >= 0.01:
        #             break
        #     if cur_it / total_it >= 0.01:
        #         break
        # if cur_it / total_it >= 0.01:
        #     break

    best_param_dict = {k: v for k, v in enumerate(best_param)}

    return [best_param_dict, min_loss]


# Below are fitted and the best params that produces the smallest loss
best_params, params_loss = fit_params(init_dict, g_exp, 5)


# simulate the best curve below
x0 = [best_params[D], 0, 0, best_params[G_b], best_params[I_b]]  # initial condition
params = list(best_params.values())


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

    dsdt = -params[0] * s
    djdt = params[0] * s - params[1] * j - params[2] * j
    J[round(t)] = j
    dldt = params[2] * time_delay(params, t, J) - params[4] * l
    dgdt = -(params[5] + params[6] * i) * g + Gprod(params, g) + params[7] * (params[1] * j + params[4] * l)
    num = params[12] ** params[11] + 1
    G_tilda = g + params[9] * (params[1] * j + params[4] * l)
    denom = (params[12] ** params[11]) * (params[14] / G_tilda) ** params[11] + 1
    didt = params[10] * params[15] * (num / denom - i / params[15])
    return [dsdt, djdt, dldt, dgdt, didt]


def time_delay(params: List[float], t: float, J: List[float]) -> float:
    if t < params[3]:
        return 0
    else:
        return J[round(t - params[3])]


def Gprod(params: List[float], g: float) -> float:
    """
    Calculate the value of G_prod at some point glucose level = g.
    """
    return params[8] / (params[8] / get_Gprod0(params) + (g - params[14]))


t = np.linspace(0, len(g_exp)-2, len(g_exp)-1)
x = odeint(gi, x0, t)

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
