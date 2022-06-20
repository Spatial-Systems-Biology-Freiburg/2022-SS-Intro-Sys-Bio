CC:=latexmk
SRCDIR:=$(CURDIR)/src/
TARGET:=$(SRCDIR)/main.tex
TARGET_HANDOUT:=handout.tex
O_DIR:=$(CURDIR)/out/
OPTIONS:=-shell-escape -pdf --output-directory=$(O_DIR) -cd

MAKEFLAGS += -j20

all: main.pdf

main.pdf:
	$(CC) $(OPTIONS) $(TARGET)

fresh: clean all

examine: examine_main

examine_main:
	pplatex -i $(O_DIR)/main.log | grep -A3 --color=always "Error\|Warning"

clean:
	rm -r $(O_DIR)/*
