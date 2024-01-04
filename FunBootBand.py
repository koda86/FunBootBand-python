import pandas as pd
import numpy as np
from scipy.linalg import svd

def initialize_variables(data, iid):
  """
  Initialize global variables
  """
  if iid is False:
    if not isinstance(data, pd.DataFrame):
      raise ValueError("Input data is not a data frame.")
      n_cluster = len(np.unique(data.columns))
      curves_per_cluster = data.shape[1] // n_cluster
      if n_cluster < 2 or n_cluster == data.shape[1]:
        raise ValueError("The header does not indicate a nested structure even though 'iid' is set to 'FALSE'.")
  else:
    # Define dummy variables (needed to avoid error in later function all)
    n_cluster = None
    curves_per_cluster = None
  
  n_time = data.shape[0]
  n_curves = data.shape[1]
  
  return n_time, n_curves, n_cluster, curves_per_cluster


# def check_arguments(band_type, alpha, iid, k_coef, B, data):
#     """
#     Function to check the validity of the input arguments
#     """
#     if not isinstance(band_type, str):
#         raise ValueError("'band_type' must be a variable of type 'character'.")
#     elif band_type not in ["confidence", "prediction"]:
#         raise ValueError("'band_type' must be either 'confidence' or 'prediction'.")
# 
#     if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
#         raise ValueError("'alpha' must be a numeric value between 0 and 1.")
# 
#     if not isinstance(iid, bool):
#         raise ValueError("'iid' must be a logical value (True or False).")
# 
#     if not isinstance(k_coef, (int, float)) or k_coef <= 0:
#         raise ValueError("'k.coef' must be a positive integer.")
# 
#     if not isinstance(B, (int, float)) or B <= 0:
#         raise ValueError("'B' must be a positive integer.")
# 
#     if data.isna().any().any():
#         raise ValueError("Function stopped due to NA's in the input data.")
#       
#     return "All checks passed successfully."


def pseudo_inverse(A, tol=np.finfo(float).eps**(2/3)):
    """
    Helper function to calculate the Moore-Penrose pseudoinverse.
    """
    U, S, Vh = svd(A)
    threshold = max(tol * S[0], 0)
    non_zero_indices = S > threshold
    S_inv = 1. / S[:len(non_zero_indices)][non_zero_indices]
    V = Vh.conj().T
    U = U[:, non_zero_indices]
    return np.dot(V[:, non_zero_indices], np.dot(np.diag(S_inv), U.conj().T))


def approximate_fourier_curves(data, k_coef, n_time, n_curves):
    """
    Function to approximate curves using Fourier functions.
    """
    time = np.arange(n_time)  # Assuming time is defined as a range from 0 to n_time-1

    # Initialize Fourier series
    fourier_s = np.ones(n_time)
    for k in range(1, k_coef * 2, 2):
        fourier_s = np.column_stack((fourier_s, np.cos(2 * np.pi * (k / 2) * time / (n_time - 1))))
        fourier_s = np.column_stack((fourier_s, np.sin(2 * np.pi * (k / 2) * time / (n_time - 1))))

    # Calculate Fourier coefficients and real values
    fourier_koeffi = np.zeros((k_coef * 2 + 1, n_curves))
    fourier_real = np.zeros((n_time, n_curves))
    for i in range(n_curves):
        fourier_koeffi[:, i] = np.dot(pseudo_inverse(np.dot(fourier_s.T, fourier_s)), np.dot(fourier_s.T, data.iloc[:, i]))
        fourier_real[:, i] = np.dot(fourier_s, fourier_koeffi[:, i])

    return fourier_koeffi, fourier_s, fourier_real


def calculate_fourier_statistics(fourier_koeffi, fourier_s, k_coef, n_curves, n_time):
    """
    Function to calculate mean Fourier curve and standard deviation.
    """
    fourier_mean = np.mean(fourier_koeffi, axis=1)
    fourier_real_mw = np.dot(fourier_s, fourier_mean)

    fourier_std1 = np.zeros((k_coef * 2 + 1, k_coef * 2 + 1, n_curves))
    for i in range(n_curves):
        diff = fourier_koeffi[:, i] - fourier_mean
        fourier_std1[:, :, i] = np.outer(diff, diff)

    fourier_cov = np.mean(fourier_std1, axis=2)
    fourier_std_all = np.sqrt(np.dot(fourier_s, np.dot(fourier_cov, fourier_s.T)))
    fourier_std = np.diag(fourier_std_all)

    return fourier_real_mw, fourier_std


def bootstrap(fourier_koeffi, fourier_s, k_coef, B, iid, n_cluster, curves_per_cluster, n_curves, n_time):
    """
    Function to perform bootstrapping on the data.
    """
    # Initialize arrays
    bootstrap_sample = np.zeros((n_time, 4))
    bootstrap_mean = np.zeros((k_coef * 2 + 1, B))
    bootstrap_real_mw = np.zeros((n_time, B))
    bootstrap_zz = np.zeros((n_curves, B), dtype=int)
    bootstrap_pseudo_koeffi = np.zeros((k_coef * 2 + 1, n_curves, B))
    bootstrap_real = np.zeros((n_time, n_curves, B))
    bootstrap_std1 = np.zeros((k_coef * 2 + 1, k_coef * 2 + 1, n_curves))
    bootstrap_cov = np.zeros((k_coef * 2 + 1, k_coef * 2 + 1, B))
    bootstrap_std_all = np.zeros((n_time, n_time, B))
    bootstrap_std = np.zeros((n_time, B))

    for i in range(B):
        if not iid:
            # Two-stage (cluster) bootstrap
            for k in range(curves_per_cluster):
                # Stage 1: Sample curve clusters with replacement
                stage_1_idx = np.random.choice(range(n_cluster), n_cluster, replace=True)
                # Stage 2: Sample within stage clusters without replacement
                curves = []
                for curve_idx in stage_1_idx:
                    curve_numbers_stage_1 = np.arange(curve_idx * curves_per_cluster - curves_per_cluster + 1,
                                                      curve_idx * curves_per_cluster + 1)
                    tmp = np.random.choice(curve_numbers_stage_1, 1, replace=False)
                    while tmp[0] in curves:
                        tmp = np.random.choice(curve_numbers_stage_1, 1)
                    curves.append(tmp[0])
                bootstrap_zz[k, i] = curves[k]
                bootstrap_pseudo_koeffi[:, k, i] = fourier_koeffi[:, bootstrap_zz[k, i]]
                bootstrap_real[:, k, i] = np.dot(fourier_s, bootstrap_pseudo_koeffi[:, k, i])
        else:
            # Ordinary (naive) bootstrap
            for k in range(n_curves):
                bootstrap_zz[k, i] = np.random.choice(n_curves, 1)[0]
                bootstrap_pseudo_koeffi[:, k, i] = fourier_koeffi[:, bootstrap_zz[k, i]]
                bootstrap_real[:, k, i] = np.dot(fourier_s, bootstrap_pseudo_koeffi[:, k, i])

        # Mean bootstrap curve and standard deviation
        bootstrap_mean[:, i] = np.mean(bootstrap_pseudo_koeffi[:, :, i], axis=1)
        bootstrap_real_mw[:, i] = np.dot(fourier_s, bootstrap_mean[:, i])

        for k in range(n_curves):
            diff = bootstrap_pseudo_koeffi[:, k, i] - bootstrap_mean[:, i]
            bootstrap_std1[:, :, k] = np.outer(diff, diff)

        bootstrap_cov[:, :, i] = np.mean(bootstrap_std1, axis=2)
        # The following line reuturns a warning because of negative numbers being used in np.sqrt().
        #  The same thing happens in the R script, and is 'tolerated' there ... why is that?
        bootstrap_std_all[:, :, i] = np.sqrt(np.dot(fourier_s, np.dot(bootstrap_cov[:, :, i], fourier_s.T)))

        for k in range(n_time):
            bootstrap_std[k, i] = bootstrap_std_all[k, k, i]

    return bootstrap_real_mw, bootstrap_std


def construct_bands(bootstrap_real_mw, bootstrap_std, fourier_real, fourier_real_mw, band_type, alpha, B, n_curves):
    """
    Function to construct statistical bands from bootstrap samples.
    """
    band_mean = np.mean(bootstrap_real_mw, axis=1)
    band_sd = np.mean(bootstrap_std, axis=1)

    if band_type == "prediction":
        cp_data = np.zeros((n_curves, B))
        cp_data_i = np.zeros((n_curves, B))

        cp_mean = 0
        cp_bound = 0
        while cp_mean < (1 - alpha):
            for i in range(B):
                for k in range(n_curves):
                    cp_data[k, i] = np.max(np.abs(fourier_real[:, k] - bootstrap_real_mw[:, i]) / bootstrap_std[:, i])
                    cp_data_i[k, i] = cp_data[k, i] < cp_bound
            cp_mean = np.mean(cp_data_i)
            cp_bound += 0.05
        band_boot = np.vstack((band_mean + cp_bound * band_sd, band_mean, band_mean - cp_bound * band_sd))

    elif band_type == "confidence":
        cc_data = np.zeros((n_curves, B))

        for i in range(B):
            for k in range(n_curves):
                cc_data[k, i] = np.max(np.abs(bootstrap_real_mw[:, i] - fourier_real_mw) / bootstrap_std[:, i])
        cc = np.quantile(cc_data, 1 - alpha)
        band_boot = np.vstack((band_mean + cc * band_sd, band_mean, band_mean - cc * band_sd))

    return band_boot


def band(data, B, alpha, iid, band_type, k_coef):
    """
    Main function to calculate statistical bands.
    """
    n_time, n_curves, n_cluster, curves_per_cluster = initialize_variables(data, iid)
    
    # Helper function to check if all arguments work. TODO: Implement try-except
    # check_arguments(band_type, alpha, iid, k_coef, B, data)
    
    # fourier_real, fourier_std = approximate_curves(data, k_coef, iid) # TODO: delete
    fourier_koeffi, fourier_s, fourier_real = approximate_fourier_curves(data, k_coef, n_time, n_curves)
    fourier_real_mw, fourier_std = calculate_fourier_statistics(fourier_koeffi, fourier_s, k_coef, n_curves, n_time)
    
    bootstrap_real_mw, bootstrap_std = bootstrap(fourier_koeffi, fourier_s, k_coef, B, iid, n_cluster, curves_per_cluster, n_curves, n_time)
    
    band_boot = construct_bands(bootstrap_real_mw, bootstrap_std, fourier_real, fourier_real_mw, band_type, alpha, B, n_curves)

    return band_boot
