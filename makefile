#MPICXX = mpicxx
#CXXFLAGS = -g -O3 -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG

MPICXX ?= mpiicpc
CXXFLAGS ?= -ipo -O3 -no-prec-div -xHost -static-intel \
	   -fopenmp \
	   -g -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG \
	   -I/usr/include/x86_64-linux-gnu \
	   -I/usr/include/i386-linux-gnu


include Makefile
