CC:=latexmk
TARGET:=main.tex
TARGET_HANDOUT:=handout.tex
OPTIONS:=-shell-escape -pdf
O_DIR:=$(CURDIR)/out/

MAKEFLAGS += -j20

all: main.pdf

main.pdf:
	$(CC) $(OPTIONS) --output-directory=$(O_DIR) $(TARGET)

fresh: clean all

examine: examine_main

examine_main:
	pplatex -i $(O_DIR)/main.log | grep -A3 --color=always "Error\|Warning"

clean:
	rm -r $(O_DIR)/*
