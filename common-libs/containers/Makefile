#CC = gcc
#CC = /opt/intel/bin/icc
ifeq ($(compiler), intel)
        CC = /opt/intel/bin/icc
else
        CC = gcc
        compiler=gcc
endif

CFLAGS = -std=c99 -O3

CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = -lcprops


all: list.o array_list.o linked_list.o
test: list.o array_list.o linked_list.o test-utils

compile: 
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c linked_list.c $(INCLUDES) $(LIBS)
	cd $(LIBS_ROOT)/commons && make

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)

linked_list.o: linked_list.h linked_list.c
	$(CC) $(CFLAGS) -c linked_list.c $(INCLUDES) $(LIBS)

test-utils:
	cd test && make compiler=$(compiler)

clean:
	rm -f *.o
