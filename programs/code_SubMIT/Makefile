# Build the MCMC inversion (finv) and the forward-model/plotting front end (ffwd) from one set of
# sources: forward_model.c holds the data, the Green's functions and the forward model for both.
#
#   make           both programs          make finv / make ffwd   one of them
#   make MPICC=mpicc MPIFC=mpif90         force the GNU wrappers (Intel oneAPI is preferred, faster)
#
# Flags. The C forward code uses fast math (v2 used -Ofast, which icx deprecates) without its
# finite-math assumption: -fno-finite-math-only keeps NaN and Inf defined, so the compiler cannot
# drop the test that rejects a degenerate proposal (cohort report 2026-10-23, item 4; PITFALLS B3).
# forward_model.c and finv.f90 also test the misfit on its IEEE bits.
MPICC ?= $(shell command -v mpiicx > /dev/null 2>&1 && echo mpiicx || \
          { command -v mpiicc > /dev/null 2>&1 && echo mpiicc || echo mpicc; })
MPIFC ?= $(shell command -v mpiifx > /dev/null 2>&1 && echo mpiifx || \
          { command -v mpiifort > /dev/null 2>&1 && echo mpiifort || \
          { command -v mpifort > /dev/null 2>&1 && echo mpifort || echo mpif90; }; })
SACHOME ?= $(shell command -v sac > /dev/null 2>&1 && dirname $$(dirname $$(command -v sac)) )
ifeq ($(strip $(SACHOME)),)
$(error SAC not found: put sac on PATH or run make SACHOME=/path/to/sac)
endif
c_optimization = -O3 -ffast-math -fno-finite-math-only
fortran_optimization = -O2
CC = $(MPICC) -mcmodel=medium $(c_optimization)
FC = $(MPIFC) -mcmodel=medium $(fortran_optimization)
libraries = -L$(SACHOME)/lib -lsac -lsacio -lm -lrt
shared_objects = forward_model.o mtdcmp.o sacsubc.o
$(info SubMIT build: MPICC=$(MPICC) MPIFC=$(MPIFC) SACHOME=$(SACHOME))

all: finv ffwd
finv: finv.o $(shared_objects)
	$(FC) finv.o $(shared_objects) $(libraries) -o finv
ffwd: ffwd.o $(shared_objects)
	$(FC) ffwd.o $(shared_objects) $(libraries) -o ffwd
finv.o: finv.f90
	$(FC) -c finv.f90
ffwd.o: ffwd.f90
	$(FC) -c ffwd.f90
mtdcmp.o: mtdcmp.f
	$(FC) -c mtdcmp.f
forward_model.o: forward_model.c sacsubc.h
	$(CC) -c forward_model.c
sacsubc.o: sacsubc.c sacsubc.h
	$(CC) -c sacsubc.c
clean:
	rm -f *.o finv ffwd
.PHONY: all clean
