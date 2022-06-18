#!/usr/bin/env python3

from dataclasses import dataclass
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


k1 = 10.0
k2 = 0.1
k3 = 4.93e-07
k4 = 80.0


def ODE(X, t):
    (A, B) = X
    return (k1 - k2*A + k3*A**2*B, k4 - k3*A**2*B)


class PDE_Solver:
    def __init__(self, dx: float, dy: float, initial_values: np.array, boundary_type: str, boundary_values, diffusion_constants: np.array, kinetics: callable):
        self.dx = dx
        self.dy = dy
        self.__values = initial_values
        
        # We will later use this matrix when calcualting the kinetics part
        self.__reducer = np.ones(initial_values.shape)
        self.__reducer[:,0,:] = 0.0
        self.__reducer[:,-1,:] = 0.0
        self.__reducer[:,:,0] = 0.0
        self.__reducer[:,:,-1] = 0.0

        self.__values_shape = initial_values.shape
        self.boundary_type = boundary_type
        self.boundary_values = (np.ones(self.__values_shape) - self.__reducer) * boundary_values
        self.kinetics = kinetics
        if len(initial_values.shape) != 3:
            raise ValueError("The initial values need to be of shape (n_chem, n_x, n_y)")
        if initial_values.shape[0] != diffusion_constants.shape[0]:
            raise ValueError("The number of provided diffusion constants has to coincide with the number of chemicals in the system.")
        self.diffusion_matrices_x = np.zeros((initial_values.shape[0], initial_values.shape[1], initial_values.shape[1]))
        self.diffusion_matrices_y = np.zeros((initial_values.shape[0], initial_values.shape[2], initial_values.shape[2]))

        self.r_x = diffusion_constants/self.dx**2
        self.r_y = diffusion_constants/self.dy**2

        self.diffusion_matrices_x = self.__initialize_diffusion_matrix(self.diffusion_matrices_x, self.r_x, boundary_type=boundary_type)
        self.diffusion_matrices_y = self.__initialize_diffusion_matrix(self.diffusion_matrices_y, self.r_y, boundary_type=boundary_type)

        if self.boundary_type == "dirichlet":
            self.__values = self.apply_dirichlet_boundary_conditions(self.__values, boundary_values)
        # elif self.boundary_type == "neumann":
        #     self.__values = self.apply_dirichlet_boundary_conditions(self.__values, boundary_values)

    def __initialize_diffusion_matrix(self, diffusion_matrix, r, boundary_type: str):
        for i in range(0, diffusion_matrix.shape[1]):
            for j in range(0, diffusion_matrix.shape[2]):
                if i==j-1 or i==j+1:
                    diffusion_matrix[:, i, j] = r
                if i==j:
                    if i==0 or i==diffusion_matrix.shape[1]-1:
                        if self.boundary_type == "dirichlet":
                            diffusion_matrix[:, i, j] = - 2.0*r
                        elif self.boundary_type == "neumann":
                            diffusion_matrix[:, i, j] = - 1.0*r
        return diffusion_matrix

    def apply_dirichlet_boundary_conditions(self, values, boundary_values):
        return self.__reducer*values + boundary_values*(np.ones(values.shape) - self.__reducer)

    def pde_rhs_dirichlet(self, X, t):
        return self.__reducer*(
            # Diffusion in x and y direction
            np.einsum('akc,abk->abc', X, self.diffusion_matrices_x) +
            np.einsum('abk,akc->abc', X, self.diffusion_matrices_y) +
            self.kinetics(X, t)
        )

    def pde_rhs_neumann(self, X, t):
        return (
            # Diffusion in x and y direction
            np.einsum('akc,abk->abc', X, self.diffusion_matrices_x) +
            np.einsum('abk,akc->abc', X, self.diffusion_matrices_y) +
            self.boundary_values +
            self.kinetics(X, t)
        )

    def ode_rhs_wrapper(self, Y, t, func):
        return self.pde_rhs_dirichlet(Y.reshape(self.__values_shape), t).flatten()

    def solve_pde(self, times):
        if self.boundary_type == "dirichlet":
            return odeint(self.ode_rhs_wrapper, self.__values.flatten(), times, args=(self.pde_rhs_dirichlet, )).reshape((len(times), )+self.__values_shape)
        elif self.boundary_type == "neumann":
            return odeint(self.ode_rhs_wrapper, self.__values.flatten(), times, args=(self.pde_rhs_neumann, )).reshape((len(times), )+self.__values_shape)


if __name__ == "__main__":
    initial_values = np.random.uniform(low=2.0, high= 10.0, size=(2, 30, 30))

    pde_solv = PDE_Solver(
        dx=20.0,
        dy=20.0,
        initial_values=initial_values,
        boundary_type="neumann",
        boundary_values=0.0,
        diffusion_constants=np.array([100.0, 5000.0]),
        kinetics=ODE
    )
    
    res = pde_solv.solve_pde(np.linspace(0, 3e-01))
    print(np.average(res, axis=(2, 3)))
    plt.imshow(res[5,0,:,:], cmap='viridis', interpolation='spline36')
    plt.show()