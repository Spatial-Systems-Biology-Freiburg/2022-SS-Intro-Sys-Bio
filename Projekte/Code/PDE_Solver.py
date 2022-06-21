#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import multiprocessing as mp
import itertools as it
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import time


def diffuse_neumann(X, t, r_x, r_y, boundary_values):
    dX = np.zeros(X.shape)
    # for i, j, k in it.product(*[range(X.shape[l]-1) for l in range(len(X.shape))]):
    #     dX[i,j,k] = (
    #         r_x[i] * (X[i,min(j+1,X.shape[1]),k] + X[i,max(j-1,0),k] - 2*X[i,j,k]) +
    #         r_y[i] * (X[i,j,min(k+1,X.shape[2])] + X[i,j,max(k-1,0)] - 2*X[i,j,k])
    #     )
    for i in range(X.shape[0]):
        if r_x[i] != 0 and r_y[i] != 0:
            dX[i,1:-1,1:-1] = (
                r_x[i] * (X[i,0:-2,1:-1] + X[i,2:,1:-1] - 2.0*X[i,1:-1,1:-1]) +
                r_y[i] * (X[i,1:-1,0:-2] + X[i,1:-1,2:] - 2.0*X[i,1:-1,1:-1])
            )
            dX[i,0,1:-1] = (
                r_x[i] * (X[i,1,1:-1] - 1.0*X[i,0,1:-1]) +
                r_y[i] * (X[i,0,2:] + X[i,0,0:-2] - 2.0*X[i,0,1:-1])
            )
            dX[i,-1,1:-1] = (
                r_x[i] * (X[i,-2,1:-1] - 1.0*X[i,-1,1:-1]) +
                r_y[i] * (X[i,-1,0:-2] + X[i,-1,2:] - 2.0*X[i,-1,1:-1])
            )
            dX[i,1:-1,0] = (
                r_y[i] * (X[i,1:-1,1] - 1.0*X[i,1:-1,0]) +
                r_x[i] * (X[i,0:-2,0] + X[i,2:,0] - 2.0*X[i,1:-1,0])
            )
            dX[i,1:-1,-1] = (
                r_y[i] * (X[i,1:-1,-2] - 1.0*X[i,1:-1,-1]) +
                r_x[i] * (X[i,0:-2,-1] + X[i,2:,-1] - 2.0*X[i,1:-1,-1])
            )
    # # Corners are still missing
    dX += boundary_values
    return dX


def diffuse_dirichlet(X, t, r_x, r_y, boundary_values):
    dX = np.zeros(X.shape)
    dX[:,1:-1,1:-1] = (
        r_x * (X[:,0:-2,1:-1] + X[:,2:,1:-1] - 2.0*X[:,1:-1,1:-1]) +
        r_y * (X[:,1:-1,0:-2] + X[:,1:-1,2:] - 2.0*X[:,1:-1,1:-1])
    )
    dX += boundary_values
    return dX


def pde_rhs_neumann(X, t, r_x, r_y, boundary_values, shape, start_pde_solve_time, kinetics):
    X = X.reshape(shape)
    print("[{: >8.4f}s] Solving ...".format(time.time()-start_pde_solve_time), end="\r")
    return diffuse_neumann(X, t, r_x, r_y, boundary_values).flatten() + kinetics(X, t).flatten()


def solve_pde(initial_values, times, dx, dy, boundary_type, boundary_values, diffusion_constants, kinetics):
    r_x = diffusion_constants/dx**2
    r_y = diffusion_constants/dy**2
    start_pde_solve_time = time.time()
    # if boundary_type == "dirichlet":
    #     res = odeint(self.ode_rhs_wrapper, initial_values.flatten(), times, args=(self.pde_rhs_dirichlet, )).reshape((len(times), )+initial_values.shape)
    if boundary_type == "neumann":
        # res = odeint(self.ode_rhs_wrapper, self.__values.flatten(), times, args=(self.pde_rhs_neumann, )).reshape((len(times), )+self.__values_shape)
        res = odeint(pde_rhs_neumann, initial_values.flatten(), times, args=(r_x, r_y, boundary_values, initial_values.shape, start_pde_solve_time, kinetics, )).reshape((len(times),) + initial_values.shape)
    print("[{: >8.4f}s] Solving Done".format(time.time()-start_pde_solve_time))
    return res


class PDE_Solver:
    def __init__(self,
        dx: float,
        dy: float,
        initial_values: np.array,
        boundary_type: str,
        boundary_values,
        diffusion_constants: np.array,
        kinetics: callable
    ):
        # Used for nice display when solving
        self.start_time = time.time()

        # Variables for the cartesian Mesh domain
        self.dx = dx
        self.dy = dy
        self.__values = initial_values
        
        # We will later use this matrix to omit certain other terms
        self.__reducer = np.ones(initial_values.shape)
        self.__reducer[:,0,:] = 0.0
        self.__reducer[:,-1,:] = 0.0
        self.__reducer[:,:,0] = 0.0
        self.__reducer[:,:,-1] = 0.0
        self.__reducer_x = self.__reducer
        self.__reducer_y = self.__reducer
        self.__reducer_x[:,:,0] = 0.0
        self.__reducer_x[:,:,-1] = 0.0
        self.__reducer_y[:,0,:] = 0.0
        self.__reducer_y[:,-1,:] = 0.0

        self.__values_shape = initial_values.shape
        self.boundary_type = boundary_type
        self.boundary_values = (np.ones(self.__values_shape) - self.__reducer) * boundary_values
        self.kinetics = kinetics
        self.dX = np.zeros(initial_values.shape)
        if len(initial_values.shape) != 3:
            raise ValueError("The initial values need to be of shape (n_chem, n_x, n_y)")
        if initial_values.shape[0] != diffusion_constants.shape[0]:
            raise ValueError("The number of provided diffusion constants has to coincide with the number of chemicals in the system.")

        self.r_x = diffusion_constants[:,np.newaxis, np.newaxis]/self.dx**2
        self.r_y = diffusion_constants[:,np.newaxis, np.newaxis]/self.dy**2
    
    def pde_rhs_dirichlet(self, X, t):
        print("[{: >8.4f}s] Solving ...".format(time.time()-self.start_pde_solve_time), end="\r")
        return self.diffuse_dirichlet(X, t, self.r_x, self.r_y, self.boundary_values) + self.kinetics(X, t)