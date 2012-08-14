#!/bin/bash
# This is an commentary line in a makefile
# Start of the makefile
FC = gfortran
PROGRAM = Henry 
#PROGRAM = Hnew

#ifeq ($(PROGRAM), Hnew)
#  SRCS = Henry_funcs.f90 Henry_sums.f90 linpacks.f90 
#  OBJS = Henry_funcs.o Henry_sums.o linpacks.o
#else
  SRCS = Henry_funcs.f90 Henry_sums.f90 LU.f90 
  OBJS = Henry_funcs.o Henry_sums.o LU.o
#endif

FCFLAGS =
OBJFLAGS = -c

all: $(PROGRAM)

$(PROGRAM): $(OBJS) 
	$(FC) $(FCFLAGS) -o $@ $@.f90 $^ 
$(OBJS): $(SRCS)
	$(FC) $(OBJFLAGS) -o $@ $*.f90
clean:
	rm *.mod *.o *.txt
# End of the makefile

