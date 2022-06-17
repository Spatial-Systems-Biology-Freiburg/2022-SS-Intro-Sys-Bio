#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


class Petri_Dish:
    def __init__(self, total_volume: float, bacteria_volume: float, n_bacteria: float, growth_rate: float, death_rate: float, initial_food: float):
        self.__current_food = initial_food
        self.__n_bacteria_current = n_bacteria
        self.__n_max = total_volume/bacteria_volume
        self.__growth_rate = growth_rate
        self.__death_rate = death_rate
    
    def advance(self, dt: float=1.0):
        dbacteria_growth = self.__growth_rate * self.__n_bacteria_current * (1 - self.__n_bacteria_current/self.__n_max)
        dbacteria_growth *= (self.__current_food/(self.__n_bacteria_current+1))
        dbacteria_growth *= (np.random.uniform(0, 1.0, 1))[0]
        dbacteria_death = self.__death_rate * self.__n_bacteria_current*np.random.uniform(0, 1.0, 1)[0]

        self.__current_food = max(self.__current_food - dt*dbacteria_growth, 0.0)
        self.__n_bacteria_current = max(self.__n_bacteria_current + dt* (dbacteria_growth - dbacteria_death), 0.0)
    
    def provide_food(self, food: float):
        if self.__current_food + food < 0.0:
            raise ValueError("Cannot provide negative food such that negative result occurs.")
        self.__current_food += food
    
    def how_much_current_food(self):
        return self.__current_food
    
    def how_much_bacteria(self):
        return self.__n_bacteria_current
    
    def take_bacteria_sample(self, n_bacteria: float):
        if n_bacteria < 0:
            raise ValueError("Taking a bacteria sample has to be a positive value.")
        self.__n_bacteria_current -= n_bacteria

    def add_bacteria(self, n_bacteria):
        if n_bacteria < 0:
            raise ValueError("Providing bacteria to Dish has to be a positive value.")
        self.__n_bacteria_current -= n_bacteria


if __name__ == "__main__":
    disch_volume = 1.0E8# The volume of the petri-dish in µm^3
    bacteria_volume = 2.0E-6#The volume of the bacteria in µm^3
    init_n_bacteria = 1E6# Initial bacteria in the Dish: no units
    init_food = 1.0E8# Initial food in the dish: in
    growth_rate = 0.1# How fast are cells reproducing (also how fast can they eat the food): in 1/s
    death_rate = 1E-1# How much % of cells are dying iteratively: in 1/s

    pd = Petri_Dish(disch_volume, bacteria_volume, init_n_bacteria, growth_rate, death_rate, init_food)

    bacteria = []
    food = []
    for i in range(0,80):
        pd.advance()
        bac = pd.how_much_bacteria()
        fod = pd.how_much_current_food()
        bacteria.append(bac)
        food.append(fod)
        

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