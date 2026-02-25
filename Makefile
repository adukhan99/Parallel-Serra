# Serra Build Makefile
CHPL_FLAGS = -M ../shared

.PHONY: all SerraNA SerraLINE clean

all: SerraNA SerraLINE

SerraNA:
	cd SerraNA && mason build -- $(CHPL_FLAGS)

SerraLINE:
	cd SerraLINE && mason build -- $(CHPL_FLAGS)

clean:
	cd SerraNA && mason clean
	cd SerraLINE && mason clean
