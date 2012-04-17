
CC = gcc
CFLAGS = -std=c99 -O3

ROOT = ~/appl/bioinfo-c

SRC_DIR = ..
REGION_DIR = $(ROOT)/bioinfo-data/features/region
COMMONS_DIR = $(ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lcprops -lcheck


all: string_utils.o list.o region_table.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o	

list.o:
	$(CC) $(CFLAGS) -c list.c $(INCLUDES) $(LIBS)

region_table.o:
	$(CC) $(CFLAGS) -c region_table.c $(INCLUDES) $(LIBS)

clean:
	rm list.o
	rm region_table.o

