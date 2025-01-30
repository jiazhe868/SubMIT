CC = mpiicx -mcmodel=medium -Ofast
FC = mpiifx -mcmodel=medium -Ofast
SACHOME=/home/staff/zjia/opt/sac-101.6a
all: finv.o mtdcmp.o sub_init.o sub_util.o sub_forward.o sub_output.o sacsubc.o alloc2d.o r8lib.o
	$(FC) finv.o mtdcmp.o sub_init.o sub_util.o sub_forward.o sub_output.o sacsubc.o alloc2d.o r8lib.o -L$(SACHOME)/lib -lsac -lsacio -lm -o finv
finv.o: finv.f90
	$(FC) -c finv.f90 -lm
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
	rm -f *.o finv
