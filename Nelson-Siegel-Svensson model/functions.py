import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from NSS import nss_rmse, constraints, nss_fit


def load_data():
    # Path to your CSV files
    risk_free_files = glob.glob("/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/Risk_free/*.csv")
    # Read and stack
    riskfree_df = pd.concat(
        [pd.read_csv(file) for file in risk_free_files],
        ignore_index=True
    )
    riskfree_df = riskfree_df.drop("2 Mo", axis=1)
    bonds_df = pd.read_csv("/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/bonds_with_timeTmaturity.csv")
    bonds_df = bonds_df[["complete_cusip", "maturity_length", "offering_date"]]
    riskfree_df = riskfree_df.rename(columns={"Date": "offering_date"})
    riskfree_df["offering_date"] = pd.to_datetime(riskfree_df["offering_date"]).dt.strftime("%Y-%m-%d")
    riskfree_df["offering_date"] = pd.to_datetime(riskfree_df["offering_date"])
    bonds_df["offering_date"] = pd.to_datetime(bonds_df["offering_date"])
    riskfree_df = riskfree_df[
        riskfree_df["offering_date"].isin(bonds_df["offering_date"])
    ]

    return riskfree_df, bonds_df


def NSS_curve(riskfree_df):
    results = []

    for date, group in riskfree_df.groupby("offering_date"):
        row = group.iloc[0].drop("offering_date")
        row = row.dropna()
        m = np.array(sorted(
            int(x.split()[0]) * (12 if "Yr" in x else 1)
            for x in row.index
            if "Mo" in x or "Yr" in x
        ))
        y = row.values
        mid_idx = len(m) // 2   # middle maturity
        three_quarter_idx = 3 * len(m) // 4   # about 75% of the way through
        print(row)
        x0 = np.array([
            y[-1],  # beta1
            y[0] - y[-1],  # beta2
            2 * y[mid_idx] - y[0] - y[-1],  # beta3
            2 * y[three_quarter_idx] - y[0] - y[-1],  # beta4
            0.0609,  # lambda1
            0.01  # lambda2
        ])

        opt = minimize(
            nss_rmse,
            x0,
            args=(y, m),
            constraints=constraints,
            method='SLSQP',
            options={
                'maxiter': 50000,
                'ftol': 1e-10,
                'disp': False
            }
        )
        plot_yield_curve(opt.x, m, y)
        results.append([
            date,
            *opt.x
        ])
    return results



def plot_yield_curve(opt_params, maturities, y_obs):
    y_fit = nss_fit(opt_params, maturities)

    plt.figure(figsize=(8, 5))

    # Observed data → only points (orange dots)
    plt.scatter(maturities, y_obs, color="orange", label="Observed")

    # Fitted NSS curve → smooth line
    plt.plot(maturities, y_fit, label="Fitted NSS")

    plt.xlabel("Maturity")
    plt.ylabel("Yield")
    plt.title("NSS Yield Curve Fit")
    plt.grid(True)
    plt.legend()

    plt.show()

