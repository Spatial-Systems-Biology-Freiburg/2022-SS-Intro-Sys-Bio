#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def RHS(y, t):
    return -0.5*y


if __name__ == "__main__":
    '''Löst eine ODE mit dem einfachen Euler-Verfahren'''
    t0 = 0.0
    tend = 10.0
    dt = 1.0
    y0 = 4.0
    
    t_vals = np.arange(t0, tend, dt)
    y_vals = np.array([y0]*len(t_vals))
    
    for i in range(len(y_vals)-1):
        y_vals[i+1] = y_vals[i] + dt * RHS(t_vals[i], y_vals[i])
      
    plt.plot(t_vals, y_vals, label="Lösung der Differentialgleichung für dt=" + str(dt))
    plt.legend()
    plt.savefig("Euler-Solutions.png")
    plt.show()
