#!/usr/bin/env python3

from PDE_Solver import PDE_Solver, save_result_plot, save_plots
import numpy as np


k1 = 10.0
k2 = 0.1
k3 = 4.93e-07
k4 = 80.0


def ODE(X, t):
    (A, B) = X
    return (k1 - k2*A + k3*A**2*B, k4 - k3*A**2*B)


if __name__ == "__main__":
    # initial_values = np.random.uniform(low=200.0, high= 900.0, size=(2, 10, 10))
    initial_values = np.zeros((2, 15, 15))
    initial_values += 900.0
    for i in range(1, 9):
        initial_values[:,i,i] = 950.0
        initial_values[:,9-i,i] = 950.0

    pde_solv = PDE_Solver(
        dx=40.0,
        dy=40.0,
        initial_values=initial_values,
        boundary_type="neumann",
        boundary_values=0.0,
        diffusion_constants=np.array([100.0, 5000.0]),
        kinetics=ODE
    )

    res = pde_solv.solve_pde(np.arange(0, 7200, 0.1))

    # Every n_th result is saved as a picture
    step = 500
    index = 0
    save_plots(res, index, step)