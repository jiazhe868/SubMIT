#FC=/opt/g95/bin/g95
FC=gfortran 
#PROG= ${HOME}/bin/teles
PROG= ./mtel3
SRCS=sub.bodyw3.f \
	sub.readray.f \
	sub.interpolation.f \
	sub.cfft.f \
	sub.clear.f \
	sub.cnvr.f \
	sub.cnvrsh.f \
	sub.instg.f \
	sub.prod.f \
	sub.qf.f \
	sub.radp3.f \
	sub.refl.f \
	sub.reflsh.f 
OBJS=${SRCS:%.f=%.o} sacio.o
LIBS=
FFLAGS = -C -O -g 


teleseis3: $(OBJS) teleseis3.f
	$(FC) $(FFLAGS) -o $(PROG) teleseis3.f $(OBJS) $(LIBS)
	gcc -o rayinformation main.c -lm 

clean:
	rm -f $(OBJS) *~ core
