# Script to test functionality/status of the functions in FunBootBand.py

import os
import pandas as pd
import FunBootBand

# Example usage
os.chdir('/home/daniel/FunBootBand-python/')
data = pd.read_csv('example_data.csv')

B = 5
alpha = 0.05
iid = True
type = 'prediction'
k_coef = 20

print(FunBootBand.check_arguments(type, alpha, iid, k_coef, B, data))

print(FunBootBand.main(data))

# Need to be translated to Python!
# plot(data[, 1], type = "l", ylim = c(-3, 3), ylab = "Amplitude")
# apply(data, 2, function(x) lines(x))
# apply(prediction.band, 1, function(x) lines(x, col = "red", lwd = 4))

