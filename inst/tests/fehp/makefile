#
# Makefile kindly by Wentian Li
#

CC      = gcc
CFLAGS  = -g -lm
 
FPMP_C  = fpmp.c 
FEHP_C  = fehp.c 
PFEHP_C = pfehp.c 
 
all:    fpmp fehp pfehp
 
fpmp:
	$(CC) $(FPMP_C) $(CFLAGS) -o fpmp
 
fehp:
	$(CC) $(FEHP_C) $(CFLAGS) -o fehp
 
pfehp:
	$(CC) $(PFEHP_C) $(CFLAGS) -o pfehp
 
clean:
	rm -f fpmp fehp pfehp
 