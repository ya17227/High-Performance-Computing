#
# Makefile Project
#

CC=mpiicc

CFLAGS= -xHOST -std=c99 -Wall -Ofast

EXE1=stencil
EXES= $(EXE1)

all: $(EXES)

$(EXES): % : %.c
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean all

clean:
	\rm -f $(EXES)
	\rm -f *.o
