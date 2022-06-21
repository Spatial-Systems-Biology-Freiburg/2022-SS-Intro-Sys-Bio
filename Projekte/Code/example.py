#!/usr/bin/env python3

from pde_functions import jpat, full_model
from pde_int import couplingMatrix, IJKth
from save_results import save_plots
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import time


def random_initialiser(xmax, ymax, NVar):
    return np.random.normal(loc=1, scale=0.1, size=xmax*ymax*NVar)


if __name__ == "__main__":
    xmax=10
    ymax=10
    NVar=7
    t_span = (0, 1000)
    t_num = 50
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
        lambda t, y: full_model(t, y, D, ind, k),
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