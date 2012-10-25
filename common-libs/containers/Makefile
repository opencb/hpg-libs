CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = -lcprops

MISC_FILES = $(LIBS_ROOT)/commons/log.c

all: list.o array_list.o
test: list.o array_list.o test-utils

compile:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS_DEBUG) -c array_list.c $(MISC_FILES) $(INCLUDES) $(LIBS)

test-utils:
	cd test && make

clean:
	rm -f list.o
	rm -f array_list.o
