import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import random
from scipy.integrate import odeint
from dataclasses import dataclass


@dataclass
class GillespieModel:
    t0: float
    y0: np.ndarray
    parameters: np.ndarray
    propensities: list
    stoichiometry: np.ndarray

    def simulate(self, steps: int):
        return gillespie(self.t0, self.y0, self.parameters, self.propensities, self.stoichiometry, steps)


def gillespie(t0, y0, parameters, propensities, stoichiometry, steps):
    # Initialize results
    t = np.zeros((steps))
    y = np.zeros((steps, len(y0)))
    # t = np.array([t0] * steps)
    # y = np.array([y0] * steps)
    t[0] = t0
    y[0] = y0

    # Begin loop
    for k in range(1, steps):
        # Calculate new rates and their sum
        rates = np.array([prop(y[k-1], t[k-1], parameters) for prop in propensities])
        a_sum = np.sum(rates)

        # Draw 2 random variables
        r1 = np.random.uniform()
        r2 = np.random.uniform()
        
        # Time for next transition
        t_next = t[k-1] + np.log(1/r1)/a_sum

        # Randomly draw a transition
        idx_trans = np.random.choice(np.arange(stoichiometry.shape[0]), p=rates/a_sum)
        y_next = np.array([
            a + b for a, b in zip(y[k-1], stoichiometry[idx_trans])
        ])

        # Update new values
        y[k] = y_next
        t[k] = t_next
    return t, y


def plot_multiple_gillespie_results(t_tot, y_tot, model):
    # Initialize plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Define the colors we want to use
    cmap = matplotlib.cm.get_cmap('viridis')
    colors = [cmap(i /( y_tot.shape[-1] - 1 )) for i in range(y_tot.shape[-1])]
    alpha = 0.2

    # Iterative over all trajectories and plot them
    for t, y in zip(t_tot, y_tot):
        # Plot all trajectories of all components
        N_comp = y_tot.shape[-1]
        for i in range(N_comp):
            ax.plot(t, y[:,i], alpha=alpha, color=colors[i])

    # Compare with exact solution
    # Setup the maximum t_range from minimum and maximum time
    t_range = np.linspace(np.min(t_tot), np.max(t_tot))
    
    # Define the ode in terms of the propensities
    ode = lambda y, t, k: np.sum([p(y, t, k)*np.array(s) for p, s in zip(model.propensities, model.stoichiometry)], axis=0)

    # Calculate the exact res by integrating the ode
    exact_res = odeint(ode, model.y0, t_range, args=(model.parameters,))

    # Plot the individual results
    for i in range(exact_res.shape[1]):
        ax.plot(t_range, exact_res[:,i], color=colors[i])

    # Make nice legend
    custom_lines =  [Line2D([0], [0], color=colors[i], lw=4, alpha=alpha) for i in range(N_comp)]
    custom_lines += [Line2D([0], [0], color=colors[i], lw=4) for i in range(N_comp)]
    labels = ["Component {} trajectories".format(i) for i in range(N_comp)]
    labels += ["Component {} exact".format(i) for i in range(N_comp)]
    ax.legend(custom_lines, labels)

    # Set labels for x and y axis
    ax.set_xlabel("Time")
    ax.set_ylabel("Component")
    plt.show()


def simulate_multiple_gillespie_runs(N_runs, model):
    # Collect results in arrays for magnitude and time
    y_tot = np.zeros((N_runs, steps, len(model.y0)))
    t_tot = np.zeros((N_runs, steps))

    # Get results from gillespie algorithm
    for l in range(N_runs):
        # t, y = gillespie(t0, y0, parameters, propensities, stoichiometry, steps)
        t, y = model_CF.simulate(steps)
        y_tot[l] = y
        t_tot[l] = t
    
    return t_tot, y_tot


def model_complex_formation():
    # Initial values
    t0 = 0.0
    y0 = np.array([25, 35, 5])

    # Parameters k0, k1 and kd
    parameters = np.array([2.0, 1, 3.0])
    
    # Define propensities
    propensities = [
        lambda y, t, k: k[0] * y[0] * y[1],
        lambda y, t, k: k[1] * y[2]
    ]

    # Every propensity (reaction) increases/decreases existing components
    stoichiometry = np.array([
        [-1, -1, 1],
        [1, 1, -1]
    ])

    return GillespieModel(t0, y0, parameters, propensities, stoichiometry)


def model_phytochromes(N_phyt: int):
    # Initial values
    t0 = 0.0
    # Define initial values randomly
    y0 = np.random.randint(low=100, high=100+2, size=(N_phyt,))
    # Alternatively set them fixed
    # y0 = np.arange(N_phyt)

    # Parameters k0, k1 and kd
    parameters = np.array([2, 1, 3])*1e-4
    
    # Define propensities foward (k0)
    propensities = [
        lambda y, t, k: k[0] *y[i] for i in range(N_phyt-1)
    ]
    # Define propensities backward (k1+kd)
    propensities += [
        lambda y, t, k: (k[1] + k[2]) * y[i+1] for i in range(N_phyt-1)
    ]

    # Every propensity (reaction) increases/decreases existing components
    # Define forward stoichiometries
    stoichiometry = [
        [0] * i + [-1, 1] + [0] * (N_phyt - i - 2) for i in range(0, N_phyt-1)
    ]
    # Define backward stoichiometries
    stoichiometry += [
        [0] * i + [1, -1] + [0] * (N_phyt - i - 2) for i in range(0, N_phyt-1)
    ]
    stoichiometry = np.array(stoichiometry)

    return GillespieModel(t0, y0, parameters, propensities, stoichiometry)


if __name__ == "__main__":
    # Define the model to solve
    model_CF = model_phytochromes(3)

    # Duration of solving
    steps = 200

    # How many gillespie runs should be displayed in the end?
    N_runs = 35

    # Simulate many gillespie runs
    t_tot, y_tot = simulate_multiple_gillespie_runs(N_runs, model_CF)
    
    # Plot the results
    plot_multiple_gillespie_results(t_tot, y_tot, model_CF)