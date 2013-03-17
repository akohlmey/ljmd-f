# -*- Makefile -*-
SHELL=/bin/sh
FC=gfortran
FFLAGS=-Wall -g -std=f95 -O3 -ffast-math -fomit-frame-pointer

# list of source files
SRC=ljmd.f90
VPATH=../src
OBJ=$(SRC:%.f90=%.o)

default: ../ljmd-serial.x

clean:
	rm -f *.mod *.o

# linker rule
../ljmd-serial.x: $(OBJ)
	$(FC) -o $@ $(FFLAGS) $^

# compilation pattern rule for objects
%.o: %.f90
	$(FC) -c $(FFLAGS) $<

# dependencies
ljmd.o: ljmd.f90
