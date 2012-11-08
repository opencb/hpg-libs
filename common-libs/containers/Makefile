CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = 

COMMONS_DIR = $(LIBS_ROOT)/commons
COMMONS_PATH = 4
MISC_OBJS = $(LIBS_ROOT)/commons/log.o

all: compile
test: list.o array_list.o test-utils

compile: compile-dependencies
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)

list.o: list.h list.c
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS_DEBUG) -c array_list.c $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(COMMONS_DIR) && make compiler=$(compiler)

test-utils:
	cd test && make

clean:
	rm list.o
	rm array_list.o
