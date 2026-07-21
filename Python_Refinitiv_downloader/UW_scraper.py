import requests
import eikon as ek
import os
import time
import importlib
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ------------------------------------------------------ README ---------------------------------------------------------- #
# this code allows, starting from a vector of bond's CUSIP, to download the underwriters of that bond's issuance, along
# with how much each underwriter partecipated and the rank description of the latter, LEAD UW, ..., ...
# ------------------------------------------------------------------------------------------------------------------------ #


# Load CUSIP list
csv_file_path = "your_csv_file_path"
df_cusip = pd.read_csv(csv_file_path)
cusip_codes = np.array(df_cusip["complete_cusip"])

# Set Eikon app key
app_key = "yourkey"
ek.set_app_key(app_key)

# Split CUSIPs in two halves
mid = len(cusip_codes) // 2
first_half = cusip_codes[:mid]
second_half = cusip_codes[mid:]

# Function to fetch data for a list of CUSIPs and save to CSV
def fetch_and_save(cusips, filename):
    all_underwriters = []

    for i, cusip in enumerate(cusips):
        try:
            bond_data, err = ek.get_data(
                instruments=[cusip],
                fields=[
                    "TR.Underwriter",
                    "TR.UnderwriterAmount",
                    "TR.UnderwriterRankDescription"
                ]
            )

            if err:
                print(f"[{i}] Error for {cusip}: {err}")
                continue

            if bond_data is None or bond_data.empty:
                print(f"[{i}] No data for {cusip}")
                continue

            # Filter only Joint Lead Managers
            bond_data = bond_data[bond_data["Underwriter Rank Description"] == "Joint Lead Manager"]
            

            for _, row in bond_data.iterrows():
                all_underwriters.append({
                    "CUSIP": cusip,
                    "Underwriter": row.get("Underwriter Name"),
                    "UnderwriterType": row.get("Underwriter Rank Description"),
                    "UnderwritingAmount": row.get("Underwriter Amount")
                })

            print(f"[{i}] Processed {cusip}, {len(bond_data)} underwriters found")

        except Exception as e:
            print(f"[{i}] Unexpected error for {cusip}: {e}")

        time.sleep(0.75)  # avoid hitting API limits
        print(bond_data)

     # Save one half
    
    df_half = pd.DataFrame(all_underwriters)
    df_half.to_csv(filename, index=False, na_rep="NA")
    


# Step 1: process first half
fetch_and_save(first_half, "/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/underwriters_half1.csv")
# Step 2: process second half
fetch_and_save(second_half, "/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/underwriters_half2.csv")

# Step 3: combine both halves and save final CSV
df1 = pd.read_csv("/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/underwriters_half1.csv")
df2 = pd.read_csv("/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/underwriters_half2.csv")
df_final = pd.concat([df1, df2], ignore_index=True)
df_final.to_csv("/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/underwriters_long_format.csv",
                index=False, na_rep="NA")

