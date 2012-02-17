#MPICXX ?= mpicxx
#CXXFLAGS ?= -Ofast -flto -fopenmp -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG -ftree-vectorizer-verbose=1
#CXXFLAGS2 ?= 
#XFLAGS ?= -fverbose-asm -masm=intel

MPICXX ?= mpiicpc
CXXFLAGS ?= -ipo -O3 -no-prec-div -xHost -static-intel \
	   -fopenmp \
	   -g -Wall -DMPICH_IGNORE_CXX_SEEK -DNDEBUG \
	   -I/usr/include/x86_64-linux-gnu \
	   -I/usr/include/i386-linux-gnu
CXXFLAGS2 ?= -vec-report1
XFLAGS ?= -fsource-asm -masm=intel


include Makefile
