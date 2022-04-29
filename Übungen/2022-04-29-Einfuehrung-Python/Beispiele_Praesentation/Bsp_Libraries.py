#!/usr/bin/env python3

# Importiere library für numerische operationen
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Benutze die libraries
    # Dokumentation von funktionen -> website der library
    x1 = np.arange(0, 7)
    x2 = np.linspace(0,7)
    # Berechne den cosinus
    y1 = np.cos(x1)
    y2 = np.cos(x2)

    # Plotte das Ergebnis
    plt.plot(x1, y1)
    plt.plot(x2, y2)
    plt.show()