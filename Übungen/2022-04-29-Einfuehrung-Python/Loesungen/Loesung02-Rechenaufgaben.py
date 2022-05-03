#!/usr/bin/env python3

import numpy as np

if __name__ == "__main__":
    # 02-01 Berechne 4^3
    print(4**3)
    
    
    # 02-02 Berechne Wurzel(2)
    print(2**0.5)
    
    
    # 02-03
    # Charlie hat 400 Melonen und will damit auf dem Mars überleben
    # Man braucht pro Tag 2.3 Melonen um genügend Nahrung zu bekommen
    # Nach {50,100,150} Tagen steigt der Verbrauch auf 2.6 Melonen pro Tag an.
    # Wie viele Tage überlebt Charlie?
    tage1 = (400 - 50*2.3)/2.6 + 50/2.3
    print(tage1)
    tage2 = (400 - 100*2.3)/2.6 + 100/2.3
    print(tage2)
    tage3 = (400 - 150*2.3)/2.6 + 150/2.3
    print(tage3)

    # Schönere Lösung
    grenzen = [50, 100, 150]
    tage = []

    for d in grenzen:
        tage.append((400 - d*2.3)/2.6 + 50/2.3)
    print(tage)
    
    
    # 02-04
    # Die Oberfläche von Hackfleisch ist ein ausgezeichneter Nährboden
    # Luca hat eine wunderschöne Hackfleischkugel mit 400g Veganem Hack
    # und einem Durchmesser von 10cm geformt und möchte
    # wissen, um wie viel sich die Oberfläche und damit die für Bakterien
    # relevante Fläche vergrößert, wenn er weitere 300g hinzufügt.
    
    # Masse ~ Volumen ~ r^3
    # Oberfläche ~ r^2
    m1 = 400
    m2 = 400 + 300
    d = 10
    O = 4*np.pi*(d/2)**2
    factor = (m2/m1)**(1/3)
    O_new = O*factor**2
    print(O)
    print(O_new)
    print(O_new/O)