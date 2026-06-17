import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize

from NSS import nss_rmse, constraints, nss_fit
from functions import load_data, NSS_curve, plot_yield_curve


def main ():
    riskfree_df, bonds_df = load_data()
    results = NSS_curve(riskfree_df)
    #bonds_rf_inter = get_risk_free(results, bonds_df)



if __name__ == "__main__":
    main()
