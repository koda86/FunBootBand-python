# Script to test functionality/status of the functions in FunBootBand.py
import os
import pandas as pd

# path_directory = os.getenv('FunBootBand_PYTHON_PATH')
path_directory = input("Enter the path to FunBootBand.py: ")
os.chdir(path_directory)

import FunBootBand
# import importlib
# importlib.reload(FunBootBand)

data = pd.read_csv('example_data.csv')
B = 10
alpha = 0.05
iid = True
band_type = "prediction"
k_coef = 50

band_boot = FunBootBand.band(data, B, alpha, iid, band_type, k_coef)

# Plot the bands
import matplotlib.pyplot as plt

plt.figure()

# Plot each curve
for column in data.columns:
    plt.plot(data[column], color='black', linewidth=1)

# Overlay the bands
for i in range(band_boot.shape[0]):
    plt.plot(band_boot[i, :], color='red', linewidth=2)

plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.show()
