CC      = gcc
CFLAGS  = -O4 -Iranlib -Iinclude -g -w
MATHLIB = -lm
SRC     = src/
SIM_OBJ = $(SRC)simmain.o \
	$(SRC)locprep.o $(SRC)getsize.o $(SRC)simped.o $(SRC)outped.o \
	$(SRC)nrutil.o
OPT_OBJ = $(SRC)diseq.o \
	$(SRC)locprep.o $(SRC)simped.o $(SRC)outped.o \
	$(SRC)nrutil.o
RANLIB  = ranlib/ranlib.a
RAN_SRC = src/ranlib.c src/com.c linpack/linpack.c
RAN_OBJ = ranlib.o com.o linpack.o

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o

all:    ranlib diseq sim

ranlib:	$(RAN_OBJ)
	(cd ranlib;ar rcv ranlib.a $(RAN_OBJ); ranlib ranlib.a)

$(RAN_OBJ): 
	(cd ranlib; $(CC) -I. -c $(RAN_SRC))

sim:    $(SIM_OBJ)
	$(CC) $(CFLAGS) $(SIM_OBJ) $(RANLIB) $(MATHLIB) -o sim

diseq:  $(OPT_OBJ)
	$(CC) $(CFLAGS) $(OPT_OBJ) $(RANLIB) $(MATHLIB) -o diseq

clean:
	rm -f $(SRC)*.o sim diseq

######################################################################
TAR_FILE = sim_tar
FILES  = Makefile $(SRC)*.c include/* ranlib/* samples/* doc/*.dvi \
	doc/*.ps
dist:
	rm -f $(TAR_FILE).Z
	tar cvf $(TAR_FILE) $(FILES)
	compress $(TAR_FILE)
######################################################################
