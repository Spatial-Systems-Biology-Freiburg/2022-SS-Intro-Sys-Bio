#!/usr/bin/env python3

if __name__ == "__main__":
    # 03-01
    # Nach einem Experiment hast to folgende Daten bekommen:
    a = [
            ["StateSize", "Proportional", "Differential", "Integral", "Total", "LastState"],
            [1, 1.500000E+02, 0.000000E+00, 0.000000E+00, 1.500000E+02, 3.000000E+01],
            [2, -8.000000E+01, -2.300000E+02, 0.000000E+00, -3.100000E+02, -1.600000E+01],
            [3, -2.250000E+02, -1.450000E+02, 0.000000E+00, -3.700000E+02, -4.500000E+01],
            [4, -2.000000E+02, 2.500000E+01, 0.000000E+00, -1.750000E+02, -4.000000E+01],
            [5, -1.950000E+02, 5.000000E+00, 0.000000E+00, -1.900000E+02, -3.900000E+01],
            [6, -1.750000E+02, 2.000000E+01, 0.000000E+00, -1.550000E+02, -3.500000E+01],
            [7, -1.600000E+02, 1.500000E+01, 0.000000E+00, -1.450000E+02, -3.200000E+01],
            [8, -1.550000E+02, 5.000000E+00, 0.000000E+00, -1.500000E+02, -3.100000E+01],
            [9, -1.250000E+02, 3.000000E+01, 0.000000E+00, -9.500000E+01, -2.500000E+01],
            [10, -7.000000E+01, 5.500000E+01, 0.000000E+00, -1.500000E+01, -1.400000E+01],
            [11, -7.000000E+01, 0.000000E+00, -4.350000E+00, -7.435000E+01, -1.400000E+01],
            [12, -5.500000E+01, 1.500000E+01, -4.533333E+00, -4.453333E+01, -1.100000E+01],
            [13, -4.500000E+01, 1.000000E+01, -4.683333E+00, -3.968333E+01, -9.000000E+00],
            [14, -4.500000E+01, 0.000000E+00, -4.833333E+00, -4.983333E+01, -9.000000E+00],
            [15, -5.000000E+00, 4.000000E+01, -4.850000E+00, 3.015000E+01, -1.000000E+00],
            [16, 5.000000E+00, 1.000000E+01, -4.833333E+00, 1.016667E+01, 1.000000E+00],
            [17, 1.500000E+01, 1.000000E+01, -4.783333E+00, 2.021667E+01, 3.000000E+00],
            [18, 3.000000E+01, 1.500000E+01, -4.683333E+00, 4.031667E+01, 6.000000E+00],
            [19, 4.000000E+01, 1.000000E+01, -4.550000E+00, 4.545000E+01, 8.000000E+00]
        ]
    # Für die Weiterverarbeitung wollen wir Benennung und Werte trennen.
    # Lege hierfür separate Variablen an und fülle sie mit den korrekten Werten
    
    # 03-02
    # Gebe die erste Spalte mit/ohne Benennung aus

    # Mit Benennung
    b1 = [sublist[0] for sublist in a[0:]]
    print(b1)

    # Ohne Benennung
    b2 = [sublist[0] for sublist in a[1:]]
    print(b2)

    # 03-03
    # Verkette die Ergebnisse von Spalte 1,2 und 3 in eine gesamte Liste
    c1 = [sublist[1] for sublist in a[1:]]
    c2 = [sublist[2] for sublist in a[1:]]
    c3 = [sublist[3] for sublist in a[1:]]
    c = c1 + c2 + c3
    print(c)

    # 03-04
    # Gib die gesamte Liste bis zur zeile 15 aus.
    print(a[0:15])

    # 03-05
    # was ist das größte Element der Liste/von jeder Spalte?
    maxes = [max(sublist) for sublist in a[1:]]
    print(maxes)