#!/usr/bin/env python3

from pde_functions import jpat
from pde_int import couplingMatrix, IJKth
from save_results import save_plots
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import time


def ODE(t, y, D, ind, k, start_time):
    print("[{: >8.4f}] Solving ...".format(time.time()-start_time), end="\r")
    dydt = np.zeros(y.shape)
    A = y[ind+0]
    B = y[ind+1]
    dydt[ind]   = k[0] - k[1]*A + k[2]*A**2*B   + k[3]*np.dot(D,A)
    dydt[ind+1] = k[4] - k[2]*A**2*B            + k[5]*np.dot(D,B)
    return dydt


if __name__ == "__main__":
    xmax=30
    ymax=30
    NVar=2
    t_span = (0, 5000)
    t_num = 1000
    y0 = np.random.uniform(low=600, high=800, size=xmax*ymax*NVar)
    bndcondition="zeroflux"
    celltype="quadratic"

    k = [
        10.0, 0.1, 4.83e-07, 100.0/40**2,
        80.0, 5000.0/40**2
    ]

    t_eval = np.linspace(t_span[0], t_span[1], t_num)
    D=couplingMatrix(xmax, ymax, bndcondition, celltype)
    ind = IJKth(1, np.arange(ymax), np.arange(xmax), ymax, NVar)

    start_time = time.time()
    print("[{: >8.4f}] Solving ...".format(0), end="\r")
    sol = solve_ivp(
        lambda t, y: ODE(t, y, D, ind, k, start_time),
        t_span,
        y0,
        method='Radau',
        jac_sparsity=jpat(D,ind),
        vectorized=True,
        t_eval=t_eval
    )
    res = sol.y.reshape((xmax, ymax, NVar, len(t_eval)))
    print("[{: >8.4f}s] Solving Done".format(time.time()-start_time))

    component_index = 0
    step = 5
    save_plots(res, component_index, step)