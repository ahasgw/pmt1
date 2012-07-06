MPICXX = mpicxx

#REPORTFLAGS = -ftree-vectorizer-verbose=1
CXXFLAGS = -g -Ofast -flto -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG \
	   -march=native -mfpmath=sse $(REPORTFLAGS)

#DIALECTFLAGS = -masm=intel
ASMFLAGS = -fverbose-asm $(DIALECTFLAGS)


include Makefile

%.s: %.cc
	$(MPICXX) $(CXXFLAGS) $(ASMFLAGS) -S $<
