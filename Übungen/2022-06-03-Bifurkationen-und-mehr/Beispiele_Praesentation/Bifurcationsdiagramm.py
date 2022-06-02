from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np



# Calculate the Right Hand Side of the ODE
def f(y, t, r):
    # UNVOLLSTÄNDIG 1
    # Wir wollen die folgende funktion ausrechnen:
    # f = ry-y^3
    return 


# Are you able to implement the same procedure with this function?
def g(y, t, a1, a2, n):
    (A, B) = y
    return a1/(1+A**n) - A, a2/(1+B**n) - B


# Can we generalize this function to be usable for f and g simultaneously?
def calculate_equilibrium(y0, r, t0, tmax):
    # Solve the ODE with values

    # UNVOLLSTÄNDIG 2
    # Hier fehlt noch etwas im Funktionsaufruf
    sol = odeint(f, y0, np.linspace(t0, tmax))

    # Test if the solution has converged
    succ = np.std(sol[-5:])/np.average(sol[-5:]) < 0.1
    
    # Return the following values:
    # y0 (initial value of y)
    # r the parameter in the ode
    # S which is the last value of the list of obtained results. We expect this to be an equilibrium
    # succ tells us if our numerical analysis determined if this was a equilibrium
    return y0, r, sol[-1,0], succ


if __name__ == "__main__":
    # Where do we start solving the ode?
    t0 = 0.0
    # Where do we stop solving the ode?
    tmax = 100.0
    
    # initial values to solve the ode for
    y_initials = np.linspace(-4, 7, num=3)
    # We vary the parameter r in the ode and plot with respect to this parameter
    r_initials = np.linspace(-4, 7, num=100)

    # Calculate the equilibria with the function we created
    points = np.array([calculate_equilibrium(yi, ri, t0, tmax) for yi in y_initials for ri in r_initials])
    # Filter them since not every calcaulated value might be a valid equilibrium
    points_filtered = np.array([p for p in points if p[-1]==True])
    # See how many we lost
    print(len(points_filtered), "of", len(points), "calculations have converged")

    # Finally plot everything to visualize
    plt.title("Bifurcationsdiagramm")
    plt.plot(points_filtered[:,1], points_filtered[:,2], '.', label="$rx - x^3$", c="k")
    plt.legend()
    plt.savefig("Bifurkationsplot.png")
    
    # UNVOLLSTÄNDIG 3
    # Warum sehen wir den plot nicht?