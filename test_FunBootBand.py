# Script to test functionality/status of the functions in FunBootBand.py
import os
import pandas as pd

os.chdir('/home/daniel/FunBootBand-python/')

data = pd.read_csv('example_data.csv')
B = 5
alpha = 0.05
iid = True
type = 'prediction'
k_coef = 20

import FunBootBand

"""
Check every function if it runs withour error
"""

print(FunBootBand.initialize_variables(data, iid))
print(FunBootBand.check_arguments(type, alpha, iid, k_coef, B, data))
print(FunBootBand.approximate_fourier_curves(data, k_coef, n_time, n_curves))
print(FunBootBand.calculate_fourier_statistics(fourier_koeffi, fourier_s, k_coef, n_curves, n_time))
print(FunBootBand.bootstrap(fourier_koeffi, fourier_s, k_coef, B, iid, n_cluster, curves_per_cluster, n_curves, n_time))

# Test the main function
print(FunBootBand.main(data))

# Needs to be translated to Python!
# plot(data[, 1], type = "l", ylim = c(-3, 3), ylab = "Amplitude")
# apply(data, 2, function(x) lines(x))
# apply(prediction.band, 1, function(x) lines(x, col = "red", lwd = 4))

