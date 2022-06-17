#!/usr/bin/env python3

from bacteria_growth import Petri_Dish
import matplotlib.pyplot as plt

class Controller():
    def __init__(self, K_p, K_i, K_d, setpoint):
        self.K_p = K_p
        self.K_i = K_i
        self.K_d = K_d
        self.setpoint = setpoint
        self.diff_hist = []
        self.integral_max = 100
    
    def difference(self, current_value):
        return self.setpoint - current_value

    def P_Controller(self, diff):
        return self.K_p * diff
    
    def D_Controller(self, diff, dt):
        return self.K_d * (diff - self.diff_hist[-1]) / dt

    def I_Controller(self, diff, dt):
        return self.K_i * sum(self.diff_hist) * dt

    def response(self, current_value, dt):
        if len(self.diff_hist) >= self.integral_max:
            self.diff_hist.pop(0)
        self.diff_hist.append(current_value)
        diff = self.difference(current_value)
        return self.P_Controller(diff) + self.D_Controller(diff, dt) + self.I_Controller(diff, dt)


if __name__ == "__main__":
    disch_volume = 1.0E8# The volume of the petri-dish in µm^3
    bacteria_volume = 2.0E-6#The volume of the bacteria in µm^3
    init_n_bacteria = 1E6# Initial bacteria in the Dish: no units
    init_food = 1.0E6# Initial food in the dish: in
    growth_rate = 0.1# How fast are cells reproducing (also how fast can they eat the food): in 1/s
    death_rate = 1E-1# How much % of cells are dying iteratively: in 1/s

    setpoint = 4e6
    K_p = 1.0e1
    K_i = 1.0e-4
    K_d = 1.0e-6

    pd = Petri_Dish(disch_volume, bacteria_volume, init_n_bacteria, growth_rate, death_rate, init_food)
    cont = Controller(K_p, K_i, K_d, setpoint)

    bacteria = []
    food = []

    add_food = 0.0

    for i in range(0,350):
        pd.advance()
        bac = pd.how_much_bacteria()
        fod = pd.how_much_current_food()

        # Now comes the controller part
        add_food = max(add_food + cont.response(current_value=bac, dt=1.0), 0.0)
        
        pd.provide_food(add_food)
        # else:
        #     pd.take_bacteria_sample(-response)
        print(add_food)
        
        bacteria.append(bac)
        food.append(fod)
        
    # Nur zum plotten
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(bacteria, color="blue", label="Bacteria")
    ax1.set_ylabel("#Bacteria")
    ax2.plot(food, color="green", linestyle="--", label="Food")
    ax2.set_ylabel("Food")
    ax1.set_title("Bacteria Growth and Food in the Disch")
    ax1.set_xlabel("Time")
    fig.legend()
    plt.show()