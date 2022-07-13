#!/usr/bin/env python3

# General Imports
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.stats import qmc
import multiprocessing as mp
import itertools as it

# Custom code Imports
from stability_filtering import lsa
from pde_functions import jpat, full_model, jac_full_model, MYC1_model, jac_MYC1_model
from pde_int import couplingMatrix, IJKth
from save_results import save_plots
from example import random_initialiser


def determine_if_pattern(res, component_index, xmax, ymax):
    y = res[:,:,:,-1]
    # Determine if a pattern was created
    y_min = np.min(y, axis=(0,1))
    y_max = np.max(y, axis=(0,1))
    y_std = np.std(y, axis=(0,1))
    y_avg = np.average(y, axis=(0,1))
    
    # Check which nodes in the space count as peaks
    filt = (y > 0.5 * (y_max - y_min) + y_min)[:,:,component_index]
    n_peaks = np.sum(filt[:,:])

    # Reject if we only have fluctuations
    if (y_std / (y_max - y_min))[component_index] > 1e-2:
        return False, n_peaks, "Rejected: y_std / (y_max - y_min) > 1e-2"
    # Reject if there are too many peaks
    if n_peaks >= y[:,:,component_index].size / 2.0:
        return False, n_peaks, "Rejected: n_peaks > xmax * ymax / 2"
    # Reject if number of peaks is not large enough
    if n_peaks < 4:
        return False, n_peaks, "Rejected: n_peaks < 4"
    # Reject if there are adjacent peaks in x or y direction
    if np.sum(filt[:-1]*filt[1:]) > 0:
        return False, n_peaks, "Rejected: Adjacent peaks"
    if np.sum(filt[:,:-1]*filt[:,1:]) > 0:
        return False, n_peaks, "Rejected: Adjacent peaks"
    # Otherwise assume that a pattern was successfully generated
    return True, n_peaks, "Not rejected"


def stability_wrapper(k, diff_func: callable, model: callable, jac_model: callable, t_span, xmax, ymax, NVar, method='Radau', error_logs_file="error.logs"):
    # Actually test the parameter configuration
    diffusion_D = diff_func(k)
    res = lsa(diffusion_D, k, model, jac_model, t_span, xmax, ymax, NVar, method='Radau', error_logs_file="error.logs")

    # If the test was successfull, store results
    return (res, diffusion_D, k, t_span, xmax, ymax, NVar)


def diffusion_func_myc1(k):
    return np.diag([1, 0, k[8], 0, 0])


def solving_wrapper(p, t_eval, model, bndcondition, celltype, component_index):
    (res, diffusion_D, k, t_span, xmax, ymax, NVar) = p
    if res == False:
        return False, -1, "Rejected by stability analysis", p
    D = couplingMatrix(xmax, ymax, bndcondition, celltype)
    ind = IJKth(1, np.arange(ymax), np.arange(xmax), ymax, NVar)
    # Solve the coupled ODEs
    sol = solve_ivp(
        lambda t, y: model(t, y, D, ind, k),
        t_span,
        y0,
        method='Radau',
        jac_sparsity = jpat(D,ind),
        vectorized = True,
        t_eval = t_eval
    )
    res = sol.y.reshape((xmax, ymax, NVar, len(t_eval)))
    msg = determine_if_pattern(res, component_index, xmax, ymax)
    return msg, p


if __name__ == "__main__":
    #######################################
    ### Values to initialize the system ###
    #######################################
    # these are independent to parameter changes
    
    # Domain size
    xmax=15
    ymax=15
    
    # Number of components of the simulated system
    NVar=5
    
    # Time values to solve for (span and steps)
    t_span = (0, 1000)
    t_num = 200
    t_eval = np.linspace(*t_span, t_num)

    # Initial values instatiated randomly
    y0 = random_initialiser(xmax, ymax, NVar)
    
    # Boundary condition of the laplace operator
    bndcondition="zeroflux"

    # Grid type
    celltype="quadratic"

    # Model equations
    model = MYC1_model
    jac_model = jac_MYC1_model

    ##########################
    ### Parameter sampling ###
    ##########################
    # use latin hypercube to sample the parameter space with boundaries
    N_param = 12
    N_samples = 2000
    p_low = [np.log10(0.05)] * N_param
    p_high = [np.log10(50.0)] * N_param

    param_sampler = qmc.LatinHypercube(N_param)
    p_sample = param_sampler.random(N_samples)
    p_sample = qmc.scale(p_sample, p_low, p_high)
    p_sample = 10 ** p_sample

    # k = np.array([
    #     0.5982, 0.1405, 2.1971, 1.1245, 0.2916, 2.3028, 0.3466, 1.7822, 0.3976, 9.9829,
    #     1.2590, 2.6202, 1.5731, 5.2625, 4.8758, 0.3196, 0.1465, 2.1453, 0.5396, 56.0520,
    #     0.5131, 0.8396, 7.8041, 1.3647
    # ])

    k1 = np.array([
        6.1823,9.3728,0.7493,35.9270,0.2009,0.1153,127.0200,604.9500,0.1164,38.5770,321.8200,0.6033
    ])

    k2 = np.array([
        4.66512072, 14.00252423, 5.6586065, 2.02561514, 0.85254995, 4.07824953,
        11.44022802, 0.26374912, 3.40301038, 20.48612917, 0.45987729, 0.06057917
    ])

    p_sample = np.concatenate(([k1, k2], p_sample))

    ##########################
    ### Stability Analysis ###
    ##########################
    # Store valid results in list
    p_sample_valid = []

    # Just for nice output
    start_time = time.time()
    
    # Define the parameters (diffusion_D and k) to test
    # diffusion_D = np.diag([k[1]*k[3], 0, 0, k[15]*k[16], k[18]*k[19], 0, 0])
    diff_func = lambda k: np.diag([1, 0, k[8], 0, 0])

    # Parallelize stability analysis
    p = mp.Pool()

    p_sample_valid = p.starmap(
        stability_wrapper,
        zip(
            p_sample,
            it.repeat(diffusion_func_myc1),
            it.repeat(model),
            it.repeat(jac_model),
            it.repeat(t_span),
            it.repeat(xmax),
            it.repeat(ymax),
            it.repeat(NVar)
        )
    )

    # Iterate over all parameters
    # for k in p_sample:
    #     p_sample_valid.append(stability_wrapper(k, diff_func, model, jac_model, t_span, xmax, ymax, NVar))

    #############################
    ### Solving the Equations ###
    #############################
    # Solve the coupled ODEs and test for patterns
    component_index = 1

    # Iterate over all valid parameter combinations
    # for p in p_sample_valid:

    res = p.starmap(
        solving_wrapper,
        zip(
            p_sample_valid,
            it.repeat(t_eval),
            it.repeat(model),
            it.repeat(bndcondition),
            it.repeat(celltype),
            it.repeat(component_index)
        )
    )
