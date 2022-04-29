#!/usr/bin/env python3


# Funktion ohne parameter
def gibt_schoene_sachen_aus():
    print("Schöne Sachen")


# Funktion mit input parametern
def sag_mir_quadrat(zahl):
    print(zahl**2)


# Funktion mit output
def gib_mir_zahl():
    return 2


# berechnet aus Dichte und Volumen die Masse
def masse(dichte, volumen):
    return dichte*volumen


if __name__ == "__main__":
    # Rufe eine funktion ohne parameter auf
    gibt_schoene_sachen_aus()

    # Rufe funktion mit input parametern auf
    sag_mir_quadrat(3)

    # Funktion mit output
    zahl = gib_mir_zahl()

    # Berechne die masse von 1000l wasser mit dichte 1kg/l
    rho = 1.0
    V = 1000.0
    print(masse(rho, V))