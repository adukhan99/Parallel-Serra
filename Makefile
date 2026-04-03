MODULES := SerraNA SerraLINE Extract Analysis WrLINE

.PHONY: all all-debug clean $(MODULES) $(MODULES:%=%-debug)

all: $(MODULES)
all-debug: $(MODULES:%=%-debug)

$(MODULES):
	cd $@ && chpl-shim mason build -- --fast

$(MODULES:%=%-debug):
	cd $(@:-debug=) && chpl-shim mason build

clean:
	for m in $(MODULES); do \
		(cd $$m && mason clean); \
	done
