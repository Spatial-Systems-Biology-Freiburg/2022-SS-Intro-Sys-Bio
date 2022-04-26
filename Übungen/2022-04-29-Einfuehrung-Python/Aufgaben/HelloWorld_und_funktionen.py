#!/usr/bin/env python3

# BEISPIELE FÜR FUNKTIONEN:

def FunktionZähltBis3():
    print("1")
    print("2")
    print("3")


def FunktionZähltBis3_advanced():
    for i in range(3):
        print(i)


def FunktionZähltBisN(N):
    for i in range(N):
        print(i)


def FunktionDieEineListeErstellt():
    # 1. Methode
    liste1 = [1, 2, 3]
    print(liste1)

    # 2. Methode
    liste2 = [i for i in range(10)]
    print(liste2)

    # 3. Methode
    liste3 = []
    for i in range(10):
        liste3.append(i)
    print(liste3)


def FunktionDieEinErgebnisAusgibt_Quadrat(x):
    return x**2


def schöneFunktionMitEinerBeschreibung():
    '''Diese Funktion ist ganz besonders toll. Schau wie sie glänzt und strahlt.'''
    print("Ich bin so schön")


if __name__ == "__main__":
    print("Hello World")

    # FunktionZähltBis3()

    # FunktionZähltBis3_advanced()

    # FunktionZähltBisN(10)

    # FunktionDieEineListeErstellt()

    # schöneFunktionMitEinerBeschreibung()

    # AUFGABEN:
    # 1) Definiere eine Funktion, die alle Zahlen zwischen N_low und N_high ausgibt
    #
    # 2) Definiere eine Funktion, welche die Fibonacci-Zahlen bis Index N berechnet
    #   - Siehe: https://en.wikipedia.org/wiki/Fibonacci_number
    #   2a) Berechne die Hasen-Population wie im Link oben erklärt nach 2, 4, 8, 16, 32 Monaten
    #   2b) und für Anfangshasen-Populationen von 0,1,2,3,4,5 Hasen
    #
    # 3) Definiere eine Funktion, die mit dem Heron-Verfahren die Wurzel einer positiven Zahl berechnet
    #   - Siehe: https://de.wikipedia.org/wiki/Heron-Verfahren
    #
    # 4) Definiere eine Funktion, die den Mittelwert einer Liste von Zahlen berechnet
