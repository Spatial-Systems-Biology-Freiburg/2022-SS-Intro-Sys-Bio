CC:=latexmk
SRCDIR:=$(CURDIR)/src/
TARGET:=main.tex
TARGET_HANDOUT:=handout.tex
O_DIR:=$(CURDIR)/out/
OPTIONS:=-pdflatex=lualatex -pdf -shell-escape --output-directory=$(O_DIR) -cd


MAKEFLAGS += -j20

all: main.pdf handout.pdf

main.pdf:
	TEXINPUTS=.//: $(CC) $(OPTIONS) $(SRCDIR)/$(TARGET)

handout.pdf:
	TEXINPUTS=.//: $(CC) $(OPTIONS) $(SRCDIR)/$(TARGET_HANDOUT)

fresh: clean all

examine: examine_main examine_handout

examine_main:
	pplatex -i $(O_DIR)/main.log | grep -A3 --color=always "Error\|Warning"

examine_handout:
	pplatex -i $(O_DIR)/handout.log | grep -A3 --color=always "Error\|Warning"

clean:
	rm -r $(O_DIR)/*
