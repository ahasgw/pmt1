MPICXX ?= mpicxx
CXXFLAGS ?= -g -O3 -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -UNDEBUG
XFLAGS2 ?=


PACKAGE = pmt0
TARGETS = pmt0
SOURCES = Makefile \
	  makefile \
	  cart.cc cart.hh \
	  conf.hh \
	  main.cc \
	  node.cc node.hh \
	  output.cc output.hh \
	  ptcl.hh \
	  random.hh \
	  signal.hh \
	  timer.hh \
	  vec.hh
LIBS = -lrt -lm

SRCS = $(filter %.cc,$(SOURCES))
HDRS = $(filter %.hh,$(SOURCES))
ASMS = $(SRCS:.cc=.s) $(SRCS:.cc=.ii)


.SUFFIXES: .cc .hh .ii .o .s
%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c $<
%.s: %.cc
	$(MPICXX) $(CXXFLAGS) $(XFLAGS) -S $<


.PHONY: all asm clean more-clean dist
all: $(TARGETS)
asm: $(ASMS)
clean:
	-rm $(TARGETS) *.o
more-clean: clean
	-rm $(ASMS)
dist:
	-[ -d /tmp/$(PACKAGE)/ ] && rm -r /tmp/$(PACKAGE)/
	mkdir /tmp/$(PACKAGE) && cp $(SOURCES) /tmp/$(PACKAGE) && \
	tar -chof - --directory /tmp $(PACKAGE) \
	| gzip -c > $(PACKAGE)_$(shell date -u +%Y%m%d.%H%M%S).tar.gz


pmt0: main.o cart.o node.o output.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
main.o: main.cc conf.hh node.hh signal.hh vec.hh timer.hh
cart.o: cart.cc cart.hh conf.hh node.hh ptcl.hh random.hh vec.hh timer.hh
node.o: node.cc cart.hh conf.hh node.hh vec.hh timer.hh
output.o: output.cc ptcl.hh vec.hh
