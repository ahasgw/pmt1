MPICXX ?= mpicxx
CXXFLAGS ?= -DMPICH_IGNORE_CXX_SEEK -UNDEBUG


PACKAGE = pmt0
TARGETS = pmt0
SOURCES = Makefile \
	  makefile.gnu makefile.intel makefile.pgi \
	  cart.cc cart.hh \
	  conf.hh \
	  main.cc \
	  node.cc node.hh \
	  output.cc output.hh \
	  ptcl.hh \
	  random.hh \
	  signal.hh \
	  timer.hh \
	  type.hh \
	  vec.hh
LIBS = -lrt -lm

SRCS = $(filter %.cc,$(SOURCES))
HDRS = $(filter %.hh,$(SOURCES))
OBJS = $(SRCS:.cc=.o)
ASMS = $(SRCS:.cc=.s)


.SUFFIXES: .cc .hh .o .s
%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c $<
%.s: %.cc
	$(MPICXX) $(CXXFLAGS) -S $<


.PHONY: all asm clean more-clean dist
all: $(TARGETS)
asm: $(ASMS)
clean:
	-rm $(TARGETS) *.o
clean-asm:
	-rm $(ASMS)
clean-all: clean clean-asm
dist:
	-[ -d /tmp/$(PACKAGE)/ ] && rm -r /tmp/$(PACKAGE)/
	mkdir /tmp/$(PACKAGE) && cp $(SOURCES) /tmp/$(PACKAGE) && \
	tar -chof - --directory /tmp $(PACKAGE) \
	| gzip -c > $(PACKAGE)_$(shell date -u +%Y%m%d.%H%M%S).tar.gz


pmt0: main.o cart.o node.o output.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
main.o: main.cc conf.hh type.hh vec.hh timer.hh node.hh signal.hh
cart.o: cart.cc cart.hh node.hh ptcl.hh type.hh vec.hh timer.hh conf.hh output.hh random.hh
node.o: node.cc node.hh conf.hh type.hh vec.hh timer.hh cart.hh ptcl.hh
output.o: output.cc output.hh ptcl.hh type.hh vec.hh
