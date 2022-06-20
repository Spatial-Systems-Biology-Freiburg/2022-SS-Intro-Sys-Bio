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

        if self.boundary_type == "dirichlet":
            self.__values = self.apply_dirichlet_boundary_conditions(self.__values, boundary_values)

        print("[{: >8.4f}s] Initialization Done".format(time.time()-self.start_time))

    def apply_dirichlet_boundary_conditions(self, values, boundary_values):
        return self.__reducer*values + boundary_values*(np.ones(values.shape) - self.__reducer)

    def diffuse_dirichlet(self, X, t):
        self.dX *= 0.0
        self.dX[:,1:-1,1:-1] = (
            self.r_x * (X[:,0:-2,1:-1] + X[:,2:,1:-1] - 2.0*X[:,1:-1,1:-1]) +
            self.r_y * (X[:,1:-1,0:-2] + X[:,1:-1,2:] - 2.0*X[:,1:-1,1:-1])
        )
        self.dX += self.boundary_values
        return self.dX

    def diffuse_neumann(self, X, t):
        self.dX *= 0.0
        self.dX[:,1:-1,1:-1] = (
            self.r_x * (X[:,0:-2,1:-1] + X[:,2:,1:-1] - 2.0*X[:,1:-1,1:-1]) +
            self.r_y * (X[:,1:-1,0:-2] + X[:,1:-1,2:] - 2.0*X[:,1:-1,1:-1])
        )
        self.dX[:,0,1:-1] = (
            self.r_x[:,0,:] * (X[:,1,1:-1] - 1.0*X[:,0,1:-1]) +
            self.r_y[:,:,0] * (X[:,0,2:] + X[:,0,0:-2] - 2.0*X[:,0,1:-1])
        )
        self.dX[:,-1,1:-1] = (
            self.r_x[:,-1,:] * (X[:,-2,1:-1] - 1.0*X[:,-1,1:-1]) +
            self.r_y[:,:,-1] * (X[:,-1,0:-2] + X[:,-1,2:] - 2.0*X[:,-1,1:-1])
        )
        self.dX[:,1:-1,0] = (
            self.r_y[:,:,0] * (X[:,1:-1,1] - 1.0*X[:,1:-1,0]) +
            self.r_x[:,:,0] * (X[:,0:-2,0] + X[:,2:,0] - 2.0*X[:,1:-1,0])
        )
        self.dX[:,1:-1,-1] = (
            self.r_y[:,:,-1] * (X[:,1:-1,-2] - 1.0*X[:,1:-1,-1]) +
            self.r_x[:,:,-1] * (X[:,0:-2,-1] + X[:,2:,-1] - 2.0*X[:,1:-1,-1])
        )
        self.dX += self.boundary_values
        return self.dX
    
    def pde_rhs_neumann(self, X, t):
        print("[{: >8.4f}s] Solving ...".format(time.time()-self.start_pde_solve_time), end="\r")
        return self.diffuse_neumann(X, t) + self.kinetics(X, t)
    
    def pde_rhs_dirichlet(self, X, t):
        print("[{: >8.4f}s] Solving ...".format(time.time()-self.start_pde_solve_time), end="\r")
        return self.diffuse_dirichlet(X, t) + self.kinetics(X, t)

    def ode_rhs_wrapper(self, Y, t, func):
        return func(Y.reshape(self.__values_shape), t).flatten()

    def solve_pde(self, times):
        self.start_pde_solve_time = time.time()
        if self.boundary_type == "dirichlet":
            res = odeint(self.ode_rhs_wrapper, self.__values.flatten(), times, args=(self.pde_rhs_dirichlet, )).reshape((len(times), )+self.__values_shape)
        elif self.boundary_type == "neumann":
            res = odeint(self.ode_rhs_wrapper, self.__values.flatten(), times, args=(self.pde_rhs_neumann, )).reshape((len(times), )+self.__values_shape)
        print("[{: >8.4f}s] Solving Done".format(time.time()-self.start_pde_solve_time))
        return res


def save_result_plot(i, u, index, min, max, start_time=None, output_folder=Path("./out/")):
    if start_time!=None:
        print("[{: >8.4f}s] Saving Plots ...".format(time.time()-start_time), end="\r")
    fig, ax = plt.subplots()
    im = ax.imshow(
        u,
        vmin=min[index],
        vmax=max[index],
        cmap='viridis',
        interpolation='spline36'
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    fig.savefig(output_folder / "image_{:08d}.png".format(i))
    plt.close(fig)


def save_plots(res, index, step, threads=None, output_folder=Path("./out/")):
    start_time = time.time()
    if threads==None:
        threads=os.cpu_count()
    with mp.Pool(threads) as p:
        p.starmap(save_result_plot, zip(
            range(0, res.shape[0], step),
            [res[i,index,:,:] for i in range(0, res.shape[0], step)],
            it.repeat(index),
            it.repeat(np.min(res, axis=(0, 2, 3))),
            it.repeat(np.max(res, axis=(0, 2, 3))),
            it.repeat(start_time),
            it.repeat(output_folder)
        ))
    print("[{: >8.4f}s] Saved all Plots".format(time.time()-start_time))