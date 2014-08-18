#
# Makefile for ranlib
#

SRC= src/
LINPACK= linpack/
OBJ= com.o ranlib.o linpack.o

CC = gcc

all: 
	$(CC) -c -I. $(SRC)ranlib.c $(SRC)com.c $(LINPACK)linpack.c
	ar rcv ranlib.a $(OBJ)
	ranlib ranlib.a
	rm -f $(OBJ)

