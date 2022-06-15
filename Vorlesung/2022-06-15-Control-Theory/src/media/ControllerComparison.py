#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def P_Functor(e: list, Kp: float):
    return Kp*e[-1]


def D_Functor(e: list, dt: float, Kd: float):
    if len(e) > 1:
        return (e[-1] - e[-2])/dt
    return 0.0


def I_Functor(e: list, dt: float, Ki: float):
    return sum(e)*dt*Ki


class Plant:
    def __init__(self):
        self.o2 = 1.0
        self.alpha = 3.0

    def apply_control_signal(self, u):
        t = np.random.random_sample()
        self.o2 += self.alpha * u*t

    def measure(self):
        return self.o2


if __name__ == "__main__":
    # Set iteration variables
    t = 0.0
    tmax = 0.2
    dt = 0.1
    
    # Define target values
    target_o2 = 1.5

    # Define Controller parameters
    Kp = 1.0
    Kd = 0.2
    Ki = 20.0

    pl = Plant()

    # Store differences
    e = []

    # Run the control loop
    while t < tmax:
        y = pl.measure()
        e.append(target_o2 - y)
        u = P_Functor(e, Kp) + D_Functor(e, dt, Kd) + I_Functor(e, dt, Ki)
        pl.apply_control_signal(u)

        t += dt
    
    plt.plot(e)
    plt.show()