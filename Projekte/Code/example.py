#!/usr/bin/env python3

from pde_functions import jpat
from pde_int import couplingMatrix, IJKth
from save_results import save_plots
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import time


def random_initialiser(xmax, ymax, NVar):
    return np.random.normal(loc=1, scale=0.1, size=xmax*ymax*NVar)


def full_model(t, y, D, ind, k, start_time):
    print("[{: >8.4f}] Solving ...".format(time.time()-start_time), end="\r")
    dydt = np.zeros(y.shape)
    TTG1 = y[ind]
    GL1  = y[ind+1]
    GL3  = y[ind+2]
    TRY  = y[ind+3]
    CPC  = y[ind+4]
    AC1  = y[ind+5]
    AC2  = y[ind+6]
    dydt[ind]   = k[0] - TTG1*(k[1] + k[2]*GL3) + (k[1]*k[3])*np.dot(D,TTG1)
    dydt[ind+1] = k[4] + k[5]*AC2 - GL1*(k[6] + k[7]*GL3)
    dydt[ind+2] = k[8] + (k[22]*k[10]*AC2*AC2)/(k[22]+AC2*AC2) \
                - GL3*(k[11] + k[7]*GL1 + k[13]*CPC) \
                + (k[23]*k[9]*AC1*AC1)/(k[23]+AC1*AC1) \
                - k[2]*TTG1*GL3 - k[12]*TRY*GL3
    dydt[ind+3] = k[14]*AC1*AC1 - TRY*k[15] - TRY*GL3*k[12] + k[15]*k[16]*np.dot(D,TRY)
    dydt[ind+4] = k[17]*AC2*AC2 - k[18]*CPC - k[13]*CPC*GL3 + k[18]*k[19]*np.dot(D,CPC)
    dydt[ind+5] = k[2]*GL3*TTG1 - k[20]*AC1
    dydt[ind+6] = k[7]*GL3*GL1 - k[21]*AC2
    return dydt



if __name__ == "__main__":
    xmax=15
    ymax=15
    NVar=7
    t_span = (0, 1000)
    t_num = 200
    y0 = random_initialiser(xmax, ymax, NVar)
    bndcondition="zeroflux"
    celltype="quadratic"

    k = [
        0.5982, 0.1405, 2.1971, 1.1245, 0.2916, 2.3028, 0.3466, 1.7822, 0.3976, 9.9829,
        1.2590, 2.6202, 1.5731, 5.2625, 4.8758, 0.3196, 0.1465, 2.1453, 0.5396, 56.0520,
        0.5131, 0.8396, 7.8041, 1.3647
    ]

    t_eval = np.linspace(t_span[0], t_span[1], t_num)

    D=couplingMatrix(xmax, ymax, bndcondition, celltype)
    ind = IJKth(1, np.arange(ymax), np.arange(xmax), ymax, NVar)

    start_time = time.time()
    print("[{: >8.4f}] Solving ...".format(0), end="\r")
    sol = solve_ivp(
        lambda t, y: full_model(t, y, D, ind, k, start_time),
        t_span,
        y0,
        method='Radau',
        jac_sparsity=jpat(D,ind),
        vectorized=True,
        t_eval=t_eval
    )
    res = sol.y.reshape((xmax, ymax, NVar, len(t_eval)))
    print("[{: >8.4f}s] Solving Done".format(time.time()-start_time))

    component_index = 1
    step = 1
    save_plots(res, component_index, step)