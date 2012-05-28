CC = gcc
CFLAGS = -std=c99 -O3

LIBS_ROOT = ..
BIOFORMATS_DIR = $(LIBS_ROOT)/bioformats
COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I . -I $(COMMONS_DIR) -I $(BIOFORMATS_DIR)/features/region -I $(BIOFORMATS_DIR)/gff
LIBS = -lcprops -lcheck


all: list.o region_table.o region_table_utils.o


list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

region_table.o:
	$(CC) $(CFLAGS) -c region_table.c $(INCLUDES) $(LIBS)

region_table_utils.o:
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c region_table_utils.c $(INCLUDES) $(LIBS)

clean:
	rm list.o
	rm region_table.o
	rm region_table_utils.o
