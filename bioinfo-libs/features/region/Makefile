
CC = gcc
CFLAGS = -std=c99 -O3

ROOT = ~/appl/bioinfo-c
SRC_DIR = .
COMMONS_DIR = $(ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(COMMONS_DIR)
LIBS = -lcheck -lcprops 


all: string_utils.o region.o

region.o: region.c
	$(CC) $(CFLAGS) -c region.c $(INCLUDES) $(LIBS)
	
string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

clean:
	rm region.o

