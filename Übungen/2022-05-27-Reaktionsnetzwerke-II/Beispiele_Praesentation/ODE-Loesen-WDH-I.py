from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

def f(y, t, a, b):
    return a - b*y

if __name__ == "__main__":
    tstart = 0.0
    tend = 10.0
    y0 = 0.0
    a = 10
    b = 0.5
    
    t = np.linspace(tstart, tend)
    results = odeint(f, y0, t, (a, b))

    plt.plot(t, results, label="Ergebnisse der gelösen ODE")
    plt.legend()
    plt.show()