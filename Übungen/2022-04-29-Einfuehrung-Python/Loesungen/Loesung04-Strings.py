#!/usr/bin/env python3

if __name__ == "__main__":
	# 04-01
	# Definiere zwei strings deiner Wahl und verkette sie zu einem
	str1 = "Das ist toll"
	str2 = " .. aber nur manchmal"
	str3 = str1 + str2
	print(str3)

	# 04-02
	# Nimm die beiden strings von zuvor und teile sie an einem Zeichen auf, das in beiden strings vorkommt
	str4 = str1.split(" ")
	print(str4)
	str5 = str2.split(" ")
	print(str5)

	# 04-03
	# Gibt den letzten buchstaben von jedem Wort aus.
	str6 = [s[-1] for s in str4]
	print(str6)
	str7 = [s for s in str5]
	print(str7)

	# 04-04
	woerter = "hallo, das hier sollen einfach ein paar ganz ganz komische und kleingeschriebene Wörter sein. Mal schauen, was so passiert!"
	# Ersetze den Anfangsbuchstaben von jedem Wort durch einen Großuchstaben
	woerter_list = woerter.split(" ")
	
	for i, w in enumerate(woerter_list):
		w = w[0].upper() + w[1:]
		woerter_list[i] = w
	
	woerter_final = " ".join(woerter_list)
	print(woerter_final)