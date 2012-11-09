#CC = gcc
#CC = /opt/intel/bin/icc
ifeq ($(COMPILER), intel)
        CC = icc
else
        CC = gcc
        COMPILER=gcc
endif

CFLAGS = -std=c99 -O3

CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = -lcprops

TEST = test

COMMONS_DIR = $(LIBS_ROOT)/commons
COMMONS_PATH = 4
MISC_OBJS = $(LIBS_ROOT)/commons/log.o

all: compile linked_list_test
test: list.o array_list.o linked_list.o test-utils

compile: compile-dependencies linked_list_test
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -c linked_list.c $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(COMMONS_DIR) && make COMPILER=$(compiler)

linked_list_test: linked_list.o $(TEST)/linked_list_test.h
	$(CC) $(CFLAGS_DEBUG) linked_list.o $(TEST)/linked_list_test.c -o bin/linked-list-test $(INCLUDES) $(LIBS)

list.o: list.h list.c
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)


linked_list.o: linked_list.h linked_list.c
	$(CC) $(CFLAGS) -c linked_list.c $(INCLUDES) $(LIBS)


test-utils:
	cd test && make COMPILER=$(compiler)

clean:
	rm -f *.o
