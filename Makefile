# -*- Makefile -*-
SHELL=/bin/sh
# compiler flags
FC=gfortran
OPT=-O3 -ffast-math -fomit-frame-pointer
FFLAGS=-Wall -g -std=f95 $(OPT)
PARALLEL=-fopenmp
# list of source files
SRC=ljmd.f90
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
OBJ_PARALLEL=$(SRC:src/%.f90=Obj-parallel/%.o)
############################################

default: ljmd-serial.x ljmd-parallel.x

ljmd-serial.x:
	$(MAKE) $(MFLAGS) -C Obj-serial

ljmd-parallel.x:
	$(MAKE) $(MFLAGS) -C Obj-parallel

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-parallel clean
	rm -f ljmd-serial.x ljmd-parallel.x


