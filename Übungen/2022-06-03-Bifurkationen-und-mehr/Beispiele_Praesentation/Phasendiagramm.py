from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


# Calculate the Right Hand Side of the ODE
def g(y, t, a1, a2, n):
    (A, B) = y
    return a1/(1+B**n) - A, a2/(1+A**n) - B


if __name__ == "__main__":
    a1 = 2.0
    a2 = 3.0
    n = 3.5

    a = np.linspace(0, 3, num=50)
    b = np.linspace(0, 4, num=50)

    A, B = np.meshgrid(a, b)

    U, V = g((A, B), 0.0, a1, a2, n)
    
    speed = np.sqrt(U**2 + V**2)

    plt.figure(figsize=(8, 8))
    plt.streamplot(A, B, U, V, linewidth=4*speed/speed.max(), color="k")
    plt.xlabel("Value A")
    plt.ylabel("Value B")
    plt.show()