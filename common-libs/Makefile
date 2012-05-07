CC = gcc
CFLAGS = -std=c99 -O3

LIBS_ROOT = ..
REGION_DIR = $(LIBS_ROOT)/bioformats/features/region
COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I . -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lcprops -lcheck


all: string_utils.o list.o region_table.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

region_table.o:
	$(CC) $(CFLAGS) -c region_table.c $(INCLUDES) $(LIBS)

clean:
	-rm list.o
	-rm region_table.o

