#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def RHS(y, t):
    return -0.5*y


if __name__ == "__main__":
    '''Löst eine ODE mit dem einfachen Euler-Verfahren'''
    t0 = 0.0
    tend = 10.0
    dt = 1.0
    y0 = 4.0

    t_vals = np.arange(t0, tend, dt)
    y_vals = odeint(RHS, y0, t_vals)
    
    plt.plot(t_vals, y_vals, label="Lösung der Differentialgleichung", marker="o")
    plt.legend()
    plt.show()