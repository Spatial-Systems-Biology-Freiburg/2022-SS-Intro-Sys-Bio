#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
	# 05-01
    # Füge eine geeignete library hinzu, um numerische berechnungen zu machen
    
    # Siehe statement oben


    # 05-02
    # Finde im Internet heraus, wie man eine zufallszahl mit dieser library generiert
    # und stecke diese dann in eine variable
    r = np.random.rand()
    print(r)

    # 05-03
    # Finde heraus, wie man Listen/arrays mit dieser library erstellt
    # und erstelle eine solche mit zahlen zwischen 0 und 10
    r_list = np.random.rand(20)*10.0
    print(r_list)

    # 05-04
    # Berechne mit dieser library cosinus, sinus und tangens der zuvor generierten werte
    # und speichere die ergebnisse in verschiedenen listen
    cos_results = np.cos(r_list)
    print(cos_results)
    sin_results = np.sin(r_list)
    print(sin_results)
    tan_results = np.tan(r_list)
    print(tan_results)

    # 05-05
    # Finde eine library, um die Ergebnisse zu plotten und mache das für alle 3 ergebnisse (in die gleiche Grafik)

    # Siehe statement oben
    plt.plot(cos_results, label="Cosinus Ergebnisse")
    plt.plot(sin_results, label="Sinus Ergebnisse")
    plt.plot(tan_results, label="Tangens Ergebnisse")
    plt.legend()
    

    # 05-06
    # Finde im Internet heraus:
    # Füge titel, skalenbenennung und Kurvenbenennung hinzu

    # Kurvenbenennung siehe oben
    plt.title("Trigonometrische Funktionen an zufälligen Werten")
    plt.xlabel("Zufallswerte")
    plt.ylabel("Trigonometrische Funktionen")
    plt.show()

    # 05-07
    # Wie kann ich alle kurven in verschiedenen "subplots" darstellen?
    # Inklusive titel, skalenbenennung + kurvenbenennung
    
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.plot(cos_results, label="Cosinus Ergebnisse")
    ax1.legend()
    ax2.plot(sin_results, label="Sininus Ergebnisse")
    ax2.legend()
    ax3.plot(tan_results, label="Tangens Ergebnisse")
    ax3.legend()
    ax1.set_title("Ergebnisse für trigonometrische Funktionen")
    ax1.set_xlabel("Achse 1 x-benennung")
    ax1.set_ylabel("Achse 1 y-benennung")
    # ... usw
    plt.show()


    # 05-08
    # Mache das gleiche für Ergebnisse der Liste von Aufgabe03
    # Automatisiere das für beliebige Listen dieser Struktur
    
    