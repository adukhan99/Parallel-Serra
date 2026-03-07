.PHONY: all SerraNA SerraLINE Extract Analysis WrLINE clean

all: SerraNA SerraLINE Extract Analysis WrLINE

SerraNA:
	cd SerraNA && chpl-shim mason build

SerraLINE:
	cd SerraLINE && chpl-shim mason build

Extract:
	cd Extract && chpl-shim mason build

Analysis:
	cd Analysis && chpl-shim mason build

WrLINE:
	cd WrLINE && chpl-shim mason build

clean:
	cd SerraNA && mason clean
	cd SerraLINE && mason clean
	cd Extract && mason clean
	cd Analysis && mason clean
	cd WrLINE && mason clean
