# VAR Analysis of fMRI data in Python

import os
import numpy as np
import pandas as pd
import scipy.io
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, acf, pacf
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale
from glob import glob

# Paths 
out_path = ""
ts_input = ""
regions = []
hcp_dir = ""

# Background parameters
tp = 1200  # number of time points
r = len(regions)  # number of regions

# Get subject IDs from file names
fc_subjects = [Path(f).stem.split("_")[0] for f in os.listdir(ts_input)]
sc_subjects = [str(s) for s in dat_SC["subject"]]
subjects = list(set(fc_subjects).intersection(sc_subjects))
subjects.sort()
n = len(subjects)

# Initialize TS and SC arrays
TS = np.zeros((n, tp, r))
SC = np.zeros((n, r, r))

# Extract FC and SC
for i, id in enumerate(subjects):
    ts_file = os.path.join(ts_input, f"{id}_TS.mat")
    data = scipy.io.loadmat(ts_file)
    timeseries = data["timeseries"]
    ts_list = []
    for j in range(4):
        reps = timeseries.shape[0] // 4
        ts = timeseries[j * reps:(j + 1) * reps, :r]
        if not np.isnan(ts).all():
            ts_list.append(ts)
    if len(ts_list) == 4:
        TS[i, :, :] = sum(ts_list) / 4
        idx = np.where(dat_SC["subject"] == int(id))[0][0]
        SC[i, :, :] = dat_SC["sift2volnorm"][idx, :, :]

# Subcortical region filtering
region_df = pd.DataFrame({"r_num": np.arange(r), "r_name": regions})
subcortical_keywords = ["Thalamus-Proper", "Hippocampus", "Amygdala"]
subcort_idx = region_df["r_name"].str.contains("|".join(subcortical_keywords))
subcort_r = region_df[subcort_idx]["r_name"].tolist()
subcort_indices = region_df[subcort_idx]["r_num"].tolist()
TS_subcort = TS[:, :, subcort_indices]
SC_subcort = SC[:, :, subcort_indices]

# Plotting Time Series of subject 0
x = pd.DataFrame(scale(TS_subcort[0, :, :]), columns=[
    "L-Thalamus", "L-Hippo", "L-Amyg",
    "R-Thalamus", "R-Hippo", "R-Amyg"
])
x.plot(title="Normalized Time Series")
plt.xlabel("Time")
plt.tight_layout()
plt.show()

# ADF test
for j in range(6):
    result = adfuller(x.iloc[:, j])
    print(f"ADF statistic for {x.columns[j]}: {result[0]}, p-value: {result[1]}")

# VAR Model
model = VAR(x)
result = model.fit(maxlags=10, ic='aic')
print(f"Selected lag order: {result.k_ar}")
print(result.summary())

# Save fitted AR coefficients for all subjects
ar1 = np.zeros((n, 6, 6))
ar2 = np.zeros((n, 6, 6))
for i in range(n):
    x_i = pd.DataFrame(scale(TS_subcort[i, :, :]), columns=x.columns)
    model_i = VAR(x_i).fit(maxlags=2)
    coefs = model_i.params.values[1:].reshape(2, 6, 6)  # lag1 and lag2
    ar1[i, :, :] = coefs[0]
    ar2[i, :, :] = coefs[1]
