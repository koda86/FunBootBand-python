# Script to test functionality/status of the functions in FunBootBand.py
import os
import pandas as pd

os.chdir('/home/daniel/FunBootBand-python/')

data = pd.read_csv('example_data.csv')
B = 5
alpha = 0.05
iid = True
band_type = 'prediction' # Should not be named 'type' in Python
k_coef = 20

import FunBootBand

"""
Check if functions run without error
"""

print(FunBootBand.initialize_variables(data, iid))
print(FunBootBand.check_arguments(type, alpha, iid, k_coef, B, data))
print(FunBootBand.approximate_fourier_curves(data, k_coef, n_time, n_curves))
print(FunBootBand.calculate_fourier_statistics(fourier_koeffi, fourier_s, k_coef, n_curves, n_time))
print(FunBootBand.bootstrap(fourier_koeffi, fourier_s, k_coef, B, iid, n_cluster, curves_per_cluster, n_curves, n_time))

# Test the main function
print(FunBootBand.main(data))

# Plot the bands
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))

# Plot each curve in the DataFrame
for column in data.columns[1:]:  # Skip the first dummy column
    plt.plot(data[column], color='black', linewidth=1)

# Overlay the bands
for i in range(band_boot.shape[0]):
    plt.plot(band_boot[i, :], color='red', linewidth=4)

# plt.title("Curves with Statistical Bands")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.show()
