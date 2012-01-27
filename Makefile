MPICXX ?= mpicxx
CXXFLAGS ?= -g -O3 -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -UNDEBUG


PACKAGE = pmt0
TARGETS = pmt0
SOURCES = Makefile \
	  main.cc \
	  conf.hh \
	  node.cc node.hh \
	  timer.cc timer.hh \
	  vec.hh
LIBS = -lrt -lm

.PHONY: all clean

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c $<

all: $(TARGETS)

pmt0: main.o node.o timer.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

main.o: main.cc conf.hh timer.hh vec.hh

node.o: node.cc node.hh vec.hh

timer.o: timer.cc timer.hh

clean:
	-rm $(TARGETS) *.o

dist:
	-[ -d /tmp/$(PACKAGE)/ ] && rm -r /tmp/$(PACKAGE)/
	mkdir /tmp/$(PACKAGE) && cp $(SOURCES) /tmp/$(PACKAGE) && \
	tar -chof - --directory /tmp $(PACKAGE) \
	| gzip -c > $(PACKAGE)_$(shell date -u +%Y%m%d.%H%M%S).tar.gz
