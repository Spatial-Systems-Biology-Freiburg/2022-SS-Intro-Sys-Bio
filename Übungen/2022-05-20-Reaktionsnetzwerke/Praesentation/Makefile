CC:=latexmk
TARGET:=main.tex
TARGET_HANDOUT:=handout.tex
OPTIONS:=-shell-escape -pdf
O_DIR:=$(CURDIR)/out/

MAKEFLAGS += -j20

all: main.pdf handout.pdf

main.pdf:
	$(CC) $(OPTIONS) --output-directory=$(O_DIR) $(TARGET)

handout.pdf:
	$(CC) $(OPTIONS) --output-directory=$(O_DIR) $(TARGET_HANDOUT)

fresh: clean all

examine: examine_main examine_handout

examine_main:
	pplatex -i $(O_DIR)/main.log | grep -A3 --color=always "Error\|Warning"

examine_handout:
	pplatex -i $(O_DIR)/handout.log | grep -A3 --color=always "Error\|Warning"

clean:
	rm -r $(O_DIR)/*
