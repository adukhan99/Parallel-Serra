# Serra Build Makefile
CHPL_FLAGS = -M ../shared

.PHONY: all SerraNA SerraLINE Extract Analysis WrLINE clean

all: SerraNA SerraLINE Extract Analysis WrLINE

SerraNA:
	cd SerraNA && mason build -- $(CHPL_FLAGS)

SerraLINE:
	cd SerraLINE && mason build -- $(CHPL_FLAGS)

Extract:
	cd Extract && mason build -- $(CHPL_FLAGS)

Analysis:
	cd Analysis && mason build -- $(CHPL_FLAGS)

WrLINE:
	cd WrLINE && mason build -- $(CHPL_FLAGS)

clean:
	cd SerraNA && mason clean
	cd SerraLINE && mason clean
	cd Extract && mason clean
	cd Analysis && mason clean
	cd WrLINE && mason clean
