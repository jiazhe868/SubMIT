# Auto-detect MPI compilers: prefer Intel (fastest), else generic MPI wrappers
# (mpicc/mpifort/mpif90 work too, just slower). Override: make MPICC=... MPIFC=...
MPICC ?= $(shell command -v mpiicx > /dev/null 2>&1 && echo mpiicx || \
          { command -v mpiicc > /dev/null 2>&1 && echo mpiicc || echo mpicc; })
MPIFC ?= $(shell command -v mpiifx > /dev/null 2>&1 && echo mpiifx || \
          { command -v mpiifort > /dev/null 2>&1 && echo mpiifort || \
          { command -v mpifort > /dev/null 2>&1 && echo mpifort || echo mpif90; }; })
# Auto-detect SAC from the sac binary in PATH; override with SACHOME=...
SACHOME ?= $(shell command -v sac > /dev/null 2>&1 && dirname $$(dirname $$(command -v sac)) )
CC = $(MPICC) -mcmodel=medium -Ofast
FC = $(MPIFC) -mcmodel=medium -Ofast
ifeq ($(strip $(SACHOME)),)
$(error SAC not found: put sac on PATH or run make SACHOME=/path/to/sac)
endif
$(info SubMIT build: MPICC=$(MPICC) MPIFC=$(MPIFC) SACHOME=$(SACHOME))
all: ffwd.o mtdcmp.o sub_init.o sub_util.o sub_forward.o sub_output.o sacsubc.o alloc2d.o r8lib.o
	$(FC) ffwd.o mtdcmp.o sub_init.o sub_util.o sub_forward.o sub_output.o sacsubc.o alloc2d.o r8lib.o -L$(SACHOME)/lib -lsac -lsacio -lm -lrt -o ffwd
ffwd.o: ffwd.f90
	$(FC) -c ffwd.f90 -lm
mtdcmp.o: mtdcmp.f
	$(FC) -c mtdcmp.f -lm
sub_init.o: sub_init.c sub_header.h
	$(CC) -c sub_init.c -lm
sub_util.o: sub_util.c sub_header.h
	$(CC) -c sub_util.c -lm
sub_forward.o: sub_forward.c sub_header.h
	$(CC) -c sub_forward.c -lm
sub_output.o: sub_output.c sub_header.h
	$(CC) -c sub_output.c -lm
alloc2d.o: alloc2d.c sub_header.h
	$(CC) -c alloc2d.c -lm
sacsubc.o: sacsubc.c sacsubc.h 
	$(CC) -c sacsubc.c -lm
r8lib.o: r8lib.c r8lib.h
	$(CC) -c r8lib.c -lm
clean:
	rm -f *.o ffwd
