#!/usr/bin/env python3

if __name__ == "__main__":
	# Erstelle eine Liste mit 7 Elementen
	a = [0, 1, 2, 3]

	# neue Liste mit veränderten Werten
	b = [elem*2 for elem in a]
	print(b)