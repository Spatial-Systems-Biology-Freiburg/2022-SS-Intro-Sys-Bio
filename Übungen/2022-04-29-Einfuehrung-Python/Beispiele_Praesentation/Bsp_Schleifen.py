#!/usr/bin/env python3

if __name__ == "__main__":
    # For-schleifen
    # Schleife über eine Liste
    a = ["das", "ist", 1, "Test"]
    for element in a:
        print(element)

    # Schleife über eine range von Zahlen
    for i in range(4):
        print(i)

    # Tricks
    # enumerate, um den index der liste zu bekommen
    for i, element in enumerate(a):
        print("Index:", i, "Inhalt:", element)