from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

def f(y, t, k1, k2):
    return (-2*k1*y[0]**2, 2*k1*y[0]**2 - k2*y[1]*y[2], - k2*y[1]*y[2])

if __name__ == "__main__":
    tstart = 0.0
    tend = 10.0
    y0 = (1.0, 0.0, 0.5)
    k1 = 0.3
    k2 = 0.5
    t = np.linspace(tstart, tend)
    
    results = odeint(f, y0, t, (k1, k2))

    for i in range(results.shape[1]):
        plt.plot(t, results[:,i], label="Komponente " + str(i))
    plt.legend()
    plt.show()