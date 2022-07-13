#!/usr/bin/env python3

import numpy as np
from scipy.integrate import solve_ivp
import multiprocessing as mp
import time

from pde_functions import jpat, full_model, MYC1_model, jac_full_model, jac_MYC1_model
from example import random_initialiser
from pde_int import couplingMatrix, IJKth


def fourier_modes(xmax, ymax):
    # Determine the possible modes:
    p = np.arange(0, xmax)*np.pi/xmax
    q = np.arange(0, ymax)*np.pi/ymax
    dpq = np.unique(np.concatenate((p, q)))
    return dpq


def log_exit_message(reason, diffusion_D, k, pde, jacobian, t_span, xmax, ymax, NVar, method, error_logs_file="error.logs"):
    print("Warning: A parameter set could not be analyzed. Writing full error log to {}".format(error_logs_file))
    f = open(error_logs_file, "a")
    message  = "____________________________________\n"
    message += "Parameter set could not be analyzed!\n"
    message += "k           = {}\n".format(k)
    message += "diffusion_D = diag({})\n".format([diffusion_D[i,i] for i in range(len(diffusion_D))])
    message += "t_span      = {}\n".format(t_span)
    message += "xmax        = {}\n".format(xmax)
    message += "ymax        = {}\n".format(ymax)
    message += "NVar        = {}\n".format(NVar)
    message += reason + "\n"
    f.write(message)
    f.close()


def lsa(diffusion_D, k, pde, jacobian, t_span, xmax, ymax, NVar, method='Radau', error_logs_file="error.logs"):
    """
    Inspects the supplied parameters for stability in time and spatial instability.
    
    Parameters
    ----------
    diffusion_D : diagonal matrix
        contains the diffusion constants in a diagonal matrix
    k : array-like
        The parameter which needs to be supplied to the ode and jacobian
    pde : callable
        Pde to inspect for stability
    jacobian : callable
        The derivative of the non-spatial reaction-kinetics part of the pde
    t_span : tuple of 2 floats
        Contains start and end point to determine the steady state
    xmax : int
        Maximum x-size of the solving domain
    ymax : int
        Maximum y-size of the solving domain
    NVar : int
        Number of variables/components in the system
    method : str
        Optional same as in scipy.integrate.solve_ivp
    error_logs_file : str
        Name of the file to store error messages.

    Returns
    ----------
    s : bool
        If stability analysis was successfull return True of False 
        depending if the parameter meets the criteria
    q : None
        If stability analysis encountered an error during the process, return None
    
    ----------
    """
    # Initialize the simulation in 2d with only 2x2 cells
    n_x = 2
    n_y = 2
    y0 = np.random.normal(loc=1, scale=0.1, size=n_x * n_y * NVar)
    bndcondition="zeroflux"
    celltype="quadratic"
    D = couplingMatrix(n_x, n_y, bndcondition, celltype)
    ind = IJKth(1, np.arange(n_y), np.arange(n_x), n_y, NVar)

    # Which fourier modes do we expect in our system?
    dpq = fourier_modes(xmax, ymax)
    # Deactivate diffusion
    D = np.zeros(D.shape)
    TS = 0
    # Obtain the steady steate of the supplied ode
    try:
        sol = solve_ivp(
            lambda t, y: pde(t, y, D, ind, k),
            t_span,
            y0,
            method='Radau',
            jac_sparsity=jpat(D,ind),
            t_eval=t_span,
            vectorized=True,
        )
    except:
        reason = "Solving of pde system encountered errors"
        log_exit_message(reason, diffusion_D, k, pde, jacobian, t_span, xmax, ymax, NVar, method)
        return None
    # Obtain the steady state via the last result
    # then average over spatial dimensions
    se = sol.y[:,-1].reshape(n_x, n_y, NVar)
    # Check if the steady state was the same at every spatial point
    ss = np.average(se, axis=(0,1))
    sd = np.std(se, axis=(0,1))
    # if not raise Error
    if np.max(sd/ss) >= 1e-6:
        reason = "Steady state could not be determined correctly. Spatial variation of component {} is {: 4.2e}".format(np.argmax(sd/ss), np.max(sd/ss))
        log_exit_message(reason, diffusion_D, k, pde, jacobian, t_span, xmax, ymax, NVar, method)
        return None
    
    # Fill the jacobian with values of the steady state and provided parameters
    J = jacobian(0, ss, k)

    # Find values with positive realpart
    evs = np.linalg.eigvals(J)
    pos_eig = evs[np.real(evs) > 0.0]
    # Reject if there are positive real eigenvalues
    if len(pos_eig) > 0:
        return False
    
    # Check spatial instability stability
    for mode in dpq:
        A = J - diffusion_D * mode**2
        evs_spat = np.linalg.eigvals(A)
        if len(evs_spat[np.real(evs_spat) >= 0.0]) > 0:
            return True
    return False


if __name__ == "__main__":
    # This script currently takes 0,0805s on average per run and has an import delay of about 0.55s
    # Tested with AMD Threadripper 3960X
    # Maximum size of the simulation
    xmax = 15
    ymax = 15
    
    # Number of variables
    NVar=5
    # Just used to calculate the steady state of the equations
    t_span = (0, 100)

    # Parameters that we want to test
    k_myc1 = np.array([
        6.1823,9.3728,0.7493,35.9270,0.2009,0.1153,127.0200,604.9500,0.1164,38.5770,321.8200,0.6033
    ])

    # Supply the coefficients of the diffusion part here. This varies from model to model!
    # full_model
    # diffusion_D = np.diag([k[1]*k[3], 0, 0, k[15]*k[16], k[18]*k[19], 0, 0])
    # MYC1_model
    diffusion_D = np.diag([1, 0, k_myc1[8], 0, 0])

    start_time = time.time()
    N_runs = 200
    for i in np.linspace(0, 100, N_runs):
        res = lsa(diffusion_D, k_myc1*i, MYC1_model, jac_MYC1_model, t_span, xmax, ymax, NVar)
    print("Total time:   {: 2.4f}".format(time.time()-start_time))
    print("Time per run: {: 2.4f}".format((time.time()-start_time)/N_runs))