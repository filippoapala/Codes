import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#factor loadings
def nss_loadings(lambdas, m):
    lambda1, lambda2 = lambdas
    L1 = np.ones(len(m))
    L2 = (1 - np.exp(-lambda1 * m)) / (lambda1 * m)
    L3 = L2 - np.exp(-lambda1 * m)
    L4 = (1 - np.exp(-lambda2 * m)) / (lambda2 * m) \
         - np.exp(-lambda2 * m)
    return np.column_stack((L1, L2, L3, L4))

def nss_fit(params, m):
    beta = params[:4]
    X = nss_loadings(params[4:6],m)
    return X @ beta

#objective function
def nss_rmse(params, y, m):
    fitted = nss_fit (params, m)
    return np.sqrt(np.mean((fitted - y) ** 2))

#constraints for the objective function
constraints = [
    {'type': 'ineq', 'fun': lambda x: x[0]},            # beta1 > 0
    {'type': 'ineq', 'fun': lambda x: x[0] + x[1]},    # beta1+beta2 > 0

    {'type': 'ineq', 'fun': lambda x: 0.15 - x[4]},    # lambda1 < 0.15
    {'type': 'ineq', 'fun': lambda x: x[4] - 0.03},    # lambda1 > 0.03

    {'type': 'ineq', 'fun': lambda x: 0.025 - x[5]},   # lambda2 < 0.025
    {'type': 'ineq', 'fun': lambda x: x[5] - 0.0075}   # lambda2 > 0.0075
]


#data for testing the model
# ----------------------------------------------------------
# NSS factor loadings
# ----------------------------------------------------------

def nss_loading(lambdas, m):

    lambda1, lambda2 = lambdas

    L1 = np.ones(len(m))

    L2 = (1 - np.exp(-lambda1 * m)) / (lambda1 * m)

    L3 = L2 - np.exp(-lambda1 * m)

    L4 = (1 - np.exp(-lambda2 * m)) / (lambda2 * m) \
         - np.exp(-lambda2 * m)

    return np.column_stack([L1, L2, L3, L4])


# ----------------------------------------------------------
# NSS fitted yields
# ----------------------------------------------------------

def nss_fit(params, m):

    beta = params[:4]

    X = nss_loading(params[4:6], m)

    return X @ beta


# ----------------------------------------------------------
# Objective function (RMSE)
# ----------------------------------------------------------

def nss_rmse(params, y, m):

    fitted = nss_fit(params, m)

    return np.sqrt(np.mean((y - fitted) ** 2))


# ----------------------------------------------------------
# Constraints
# ----------------------------------------------------------

constraints = [

    {'type': 'ineq', 'fun': lambda x: x[0]},            # beta1 > 0
    {'type': 'ineq', 'fun': lambda x: x[0] + x[1]},    # beta1+beta2 > 0

    {'type': 'ineq', 'fun': lambda x: 0.15 - x[4]},    # lambda1 < 0.15
    {'type': 'ineq', 'fun': lambda x: x[4] - 0.03},    # lambda1 > 0.03

    {'type': 'ineq', 'fun': lambda x: 0.025 - x[5]},   # lambda2 < 0.025
    {'type': 'ineq', 'fun': lambda x: x[5] - 0.0075}   # lambda2 > 0.0075
]

