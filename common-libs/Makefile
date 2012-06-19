CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g


LIBS_ROOT = $(PWD)/..
BIOFORMATS_DIR = $(LIBS_ROOT)/bioformats
COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I . -I $(COMMONS_DIR) -I $(BIOFORMATS_DIR)/features/region -I $(BIOFORMATS_DIR)/gff
LIBS = -lm -lcprops

MISC_FILES = $(COMMONS_DIR)/log.c

all: list.o array_list.o region_table.o region_table_utils.o


list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

array_list.o: array_list.h array_list.c
	$(CC) $(CFLAGS_DEBUG) -c array_list.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region_table.o:
	$(CC) $(CFLAGS) -c region_table.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region_table_utils.o:
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c region_table_utils.c $(INCLUDES) $(LIBS)

clean:
	rm list.o
	rm array_list.o
	rm region_table.o
	rm region_table_utils.o
