UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    CC=gcc-6 -fopenmp
endif

ifeq ($(UNAME_S),Linux)
    CC=gcc -fopenmp
endif

CFLAGS=-O3 -march=native -std=gnu99
EXEC=test

.PHONY: all clean

all:
	${CC} matrix.c test.c ${CFLAGS} -o test

clean:
	rm -f ${EXEC}

