# Makefile for EST Clustering program using GAs.

CC=g++
CFLAGS= -Wall -g -O3 $(INC_DIRS)

GA_INC_DIR= /usr/local/include
GA_LIB_DIR= /usr/local/lib

INC_DIRS= -I$(GA_INC_DIR)
LIB_DIRS= -L$(GA_LIB_DIR)

OBJ= gaest.o dynamic.o dna.o estest.o exest.o
EXEC= gaest estest exest

dna.o: dna.cpp dna.h
	$(CC) $(CFLAGS) -c -o dna.o dna.cpp

dynamic.o: dna.h dynamic.h dynamic.cpp
	$(CC) $(CFLAGS) -c -o dynamic.o dynamic.cpp

gaest.o: gaest.cpp dynamic.h dna.h
	$(CC) $(CFLAGS) -c -o gaest.o gaest.cpp

estest.o: estest.cpp dynamic.h dna.h
	$(CC) $(CFLAGS) -c -o estest.o estest.cpp

exest.o: exest.cpp dynamic.h dna.h
	$(CC) $(CFLAGS) -c -o exest.o exest.cpp

gaest: gaest.o dynamic.o dna.o
	$(CC) -o gaest gaest.o dynamic.o dna.o $(LIB_DIRS) -lga -lm

estest: dna.o dynamic.o estest.o
	$(CC) -o estest estest.o dna.o dynamic.o

exest: exest.o dna.o dynamic.o
	$(CC) -o exest exest.o dna.o dynamic.o

clean:
	rm -f $(OBJ)

clobber:
	rm -f $(OBJ) $(EXEC)
