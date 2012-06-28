CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = $(PWD)/..

INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = -lcprops

INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include
LIBS_STATIC = $(LIBS_ROOT) -lcprops

MISC_FILES = $(LIBS_ROOT)/commons/log.c

all: list.o array_list.o region_table.o region_table_utils.o

compile:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -c region_table.c $(MISC_FILES) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c region_table_utils.c $(INCLUDES) $(LIBS)

compile-static:
	$(CC) $(CFLAGS) -c array_list.c $(INCLUDES_STATIC) $(LIBS_STATIC)
	$(CC) $(CFLAGS) -c list.c $(INCLUDES_STATIC) $(LIBS_STATIC)
	$(CC) $(CFLAGS) -c region_table.c $(MISC_FILES) $(INCLUDES_STATIC) $(LIBS_STATIC)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c region_table_utils.c $(INCLUDES_STATIC) $(LIBS_STATIC)

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS_DEBUG) -c array_list.c $(MISC_FILES) $(INCLUDES)

region_table.o:
	$(CC) $(CFLAGS) -c region_table.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region_table_utils.o:
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c region_table_utils.c $(INCLUDES) $(LIBS)

clean:
	rm list.o
	rm array_list.o
	rm region_table.o
	rm region_table_utils.o
