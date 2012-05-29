CC = gcc
CFLAGS  = -pedantic -Wall -W -O2
COMP    = $(CC) $(CFLAGS)

mofext: mofext.o matrix.o mystrings.o extend.o similar.o getline.o
	$(COMP) mofext.o matrix.o mystrings.o extend.o similar.o getline.o -o mofext
mofext.o: matrix.h mofext.c
	$(COMP) mofext.c -c
matrix.o: matrix.h matrix.c
	$(COMP) matrix.c -c
similar.o: similar.h similar.c extend.c extend.h
	$(COMP) similar.c -c
extend.o: extend.h extend.c
	$(COMP) extend.c -c
mystrings.o: mystrings.h mystrings.c
	$(COMP) mystrings.c -c
getline.o: getline.h getline.c
	$(COMP) getline.c -c
clean:
	rm -f *~
	rm -f *.o
	rm -f core*
install:
	cp mofext /usr/local/bin
uninstall:
	rm -f /usr/local/bin/mofext
