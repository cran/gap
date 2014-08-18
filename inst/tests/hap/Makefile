VERSION = 0.2.1
CMP = gcc -c
LD = gcc -lm -o
UNIFORM_RANDOM = drand48
RANDOM_SEED = srand48

all: hap mia

hap: hap.o hapfun.o cline.o 
	$(LD) hap hap.o hapfun.o cline.o 

hap.o : hap.c hap.h cline.h 
	$(CMP)  -DSEED=$(RANDOM_SEED) hap.c

hapfun.o : hapfun.c  hap.h
	$(CMP) -DURAN=$(UNIFORM_RANDOM) hapfun.c

mia: mia.o cline.o
	$(LD) mia mia.o cline.o

mia.o: mia.c
	$(CMP) mia.c

cline.o : cline.c cline.h
	$(CMP) cline.c

clean :
	rm -f *.o *~ hap-$(VERSION).tar.gz hap-$(VERSION).zip

distribution:
	mkdir hap-$(VERSION)
	cp *.c *.h *.txt test.dat Makefile* 4snps.* hap-$(VERSION)/
	rm -f hap-$(VERSION).tar.gz hap-$(VERSION).zip
	tar cvf hap-$(VERSION).tar hap-$(VERSION)/*
	gzip hap-$(VERSION).tar
	zip -ll hap-$(VERSION).zip hap-$(VERSION)/*
	rm -f hap-$(VERSION)/*
	rmdir hap-$(VERSION)





