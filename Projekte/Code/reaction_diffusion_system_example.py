#!/usr/bin/env python3

from PDE_Solver import PDE_Solver, solve_pde, save_plots
import numpy as np
import matplotlib.pyplot as plt


k1 = 10.0
k2 = 0.1
k3 = 4.93e-07
k4 = 80.0


def ODE(X, t):
    (A, B) = X
    return (k1 - k2*A + k3*A**2*B, k4 - k3*A**2*B)


k = [0.5982, 0.1405, 2.1971, 1.1245, 0.2916, 2.3028, 0.3466, 1.7822, 0.3976, 9.9829,\
        1.2590, 2.6202, 1.5731, 5.2625, 4.8758, 0.3196, 0.1465, 2.1453, 0.5396, 56.0520,\
        0.5131, 0.8396, 7.8041, 1.3647]


def full_model(X, t):
	dydt = np.zeros(X.shape)
	(TTG1, GL1, GL3, TRY, CPC, AC1, AC2) = X
	dydt[0]   = k[0] - TTG1*(k[1] + k[2]*GL3)
	dydt[0+1] = k[4] + k[5]*AC2 - GL1*(k[6] + k[7]*GL3)
	dydt[0+2] = k[8] + (k[22]*k[10]*AC2*AC2)/(k[22]+AC2*AC2) \
                - GL3*(k[11] + k[7]*GL1 + k[13]*CPC) \
                + (k[23]*k[9]*AC1*AC1)/(k[23]+AC1*AC1) \
                - k[2]*TTG1*GL3 - k[12]*TRY*GL3
	dydt[0+3] = k[14]*AC1*AC1 - TRY*k[15] - TRY*GL3*k[12]
	dydt[0+4] = k[17]*AC2*AC2 - k[18]*CPC - k[13]*CPC*GL3
	dydt[0+5] = k[2]*GL3*TTG1 - k[20]*AC1
	dydt[0+6] = k[7]*GL3*GL1 - k[21]*AC2
	return dydt


if __name__ == "__main__":
    # initial_values = np.random.uniform(low=200.0, high= 900.0, size=(2, 10, 10))
    initial_values = np.zeros((7, 15, 15))
    initial_values += 900.0
    for i in range(1, 9):
        initial_values[:,i,i] = 950.0
        initial_values[:,9-i,i] = 950.0

    dx = 40.0
    dy = 40.0
    boundary_values=0.0
    diffusion_constants=np.array([
        k[1]*k[3],
        0,
        0,
        k[15]*k[16],
        k[18]*k[19],
        0,
        0
    ])
    times = np.arange(0, 7200, 0.1)

    res = solve_pde(initial_values, times, dx, dy, "neumann", boundary_values, diffusion_constants, kinetics=full_model)

    # Every n_th result is saved as a picture
    step = 500
    index = 0
    save_plots(res, index, step)

    plt.imshow(res[-1,0,:,:])
    plt.show()