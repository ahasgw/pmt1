#MPICXX ?= mpicxx
#CXXFLAGS ?= -g -O3 -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG
#XFLAGS ?= -save-temps -fverbose-asm -masm=intel

MPICXX ?= mpiicpc
CXXFLAGS ?= -ipo -O3 -no-prec-div -xHost -static-intel \
	   -fopenmp \
	   -g -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG \
	   -I/usr/include/x86_64-linux-gnu \
	   -I/usr/include/i386-linux-gnu
XFLAGS ?= -fsource-asm -masm=intel


include Makefile
