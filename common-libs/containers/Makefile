CC = gcc
#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3

CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = 

all: list.o array_list.o
test: list.o array_list.o test-utils

compile: 
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)
	cd $(LIBS_ROOT)/commons && make

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS_DEBUG) -c array_list.c $(INCLUDES) $(LIBS)

test-utils:
	cd test && make

clean:
	rm list.o
	rm array_list.o
