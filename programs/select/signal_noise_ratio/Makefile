# SNR tool used by step 1 (selectSNR.pl). Links the SAC I/O library.
SACHOME ?= $(shell command -v sac > /dev/null 2>&1 && dirname $$(dirname $$(command -v sac)))
CC ?= gcc
all: ../snr
../snr: snr.c
	$(CC) -O2 -I$(SACHOME)/include -o ../snr snr.c -L$(SACHOME)/lib -lsacio -lm
clean:
	rm -f ../snr *.o
