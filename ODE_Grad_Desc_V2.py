import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import autograd
from autograd.builtins import tuple
import autograd.numpy as np
from scipy.integrate import odeint
G_b=5.94444444
I_b=0.00000017
def f(y,t,params):
    '''Function describing dynamics of the system'''
    #params is an array with the following parameters in order: k_js, k_gj, k_jl, k_gl, k_paramsg, k_paramsgi, eta, beta, gamma, f_gj,k_paramsi,k_lambda
    #params(0)=k_js
    #params(1)=k_gj
    #params(2)=k_jl
    #params(3)=k_gl
    #params(4)=k_xg
    #params(5)=k_xgi
    #params(6)=eta
    #params(7)=beta
    #params(8)=gamma
    #params(9)=f_gj
    #params(10)=k_paramsi
    #params(11)=k_lambda
    s, j, l, g, i = y
    #basic implementation of differential equations; time delay set to 0
    ds = -params[0] * s
    dj = params[0] * s - params[1] * j - params[2] * j
    dl = params[2] * j - params[3] * l
    g_prod0 = params[5] * G_b * I_b + params[4] * G_b
    dg = -(params[4] + params[5] * i) * g + (params[11] / (params[11] / g_prod0 + (g - G_b))) + params[6] * (params[1] * j + params[3] * l)
    num = (params[7] ** params[8]) + 1
    g_tilda = g + params[9] * (params[2] * j + params[3] * l)
    denom = (params[7] ** params[8]) * ((G_b / g_tilda) ** params[8]) + 1 #problematic for some reason
    di = params[10] * I_b * (num / denom - i / I_b)
    return np.array([ds, dj, dl, dg, di])

# Jabobian matrix wrt y
J = autograd.jacobian(f, argnum=0) # Here 0 means position 0 in function f(y, t, params) -- taking Jacobian wrt f[0] = y
# Gradient wrt to the all of the paramaters
grad_f_params = autograd.jacobian(f, argnum=2)


def ODESYS(Y,t,params, varnum):

    #Y will be length 10.
    #Y[0], Y[1], Y[2], Y[3], Y[4] are the ODEs
    #Y[5], Y[6], Y[7], Y[8], Y[9] are the sensitivities

    #ODE
    dy_dt = f(Y[0:5],t,params)#dy_dt will be a 1x5 1-d array containing the change for each variable in y
    #Sensitivities
    sens=J(Y[:5], t, params)@Y[-5::] #sense will be a 1x5 1-d array with the sensitivity for each variable of y
    # because we are dealing with 12 parameters, grad_f_params(Y[:5], t, params) will give us a 5x12 matrix with each cell expressing the impact of the parameter (column) on the y variable (row
    #to add the grad_f_param value with the 1x5 sensitivity array, we will argue that the sensitivty of the 5 y variables is the same for each parameter
    #thus we will create a 2d 5x12 array with sens duplicated 12 times as column vectors
    grad_y_params = grad_f_params(Y[:5], t, params) + np.array([sens]*12).transpose()
    #concatenate to form 12x10 matrix (first 5 columns for ode and next 5 for sensitivty), rows represent paramaters
    comp_matrix=np.concatenate([np.array([dy_dt]*12),grad_y_params.transpose()],axis=1)
    #depending on requested parameter (for whichever is being actively fitted), return accordingly
    return comp_matrix[varnum,:]

def Cost(y_obs):
    def cost(Y):
        '''Squared Error Loss'''
        n = y_obs.shape[0]
        g_modeled=Y[:, 3]
        err = np.linalg.norm(y_obs - g_modeled, 2)
        return np.sum(err) / n

    return cost
#reading data to be fitted
freader=pd.read_csv('2020-07-2021_lunch_peak.csv')
gdata=freader.loc[:,"glucose"].to_numpy()
#gdata will have glucose levels from minuts 1 to 120
gdata=gdata[:120] #truncated data to 120 for simplicity

#t will be minutes from 1 to 120
t = np.linspace(1,120,120)

#Grad Descent Implementation
params=np.array([0.1, 0.1, 0.0316, 0.1, 0.01, 0.000000316, 0.01, 100., 5., 10., 0.01, 0.0316])
cost = Cost(gdata)
grad_C=autograd.grad(cost)
#for debugging and computational simplicity, each paramater will have 3 iterations to be fitted
maxiter=3
#leanring rate set at 0.1 but can increase if needed
learning_rate=0.1
#we will modify each parameter value independently
for i in range(params.size):
    print(i)
    for j in range(maxiter):
        # strange bug where at parameter 5 (index=4), iteration 3? a degenrate solution appears
        #Solve ODE using odeint
        sol = odeint(ODESYS, y0=np.array([417., 0., 0., G_b, I_b, 0.0, 0.0, 0.0, 0.0, 0.0]), t=t, args=tuple([params,i]))
        #isolates the values of y variables
        Y = sol[:,:5]
        #cost function tracker, should observe a strictly decreasing behavior (and we do!)
        current_error = np.linalg.norm(Y[:,3] - gdata, 2)
        print(current_error)
        #parameter in question is changed as needed
        params[i]-=learning_rate*(grad_C(Y)*sol[:,-5:]).sum()
        #graph newly fitted plot, technically one iteration old
        t = np.linspace(1, 120, 120)
        plt.scatter(t, gdata, marker='.', alpha=0.5, label='G')
        plt.scatter(t, sol[:, 3], marker='.', alpha=0.5, label='estimate')
        plt.legend()
        plt.show()

print(params)


