
CC = gcc
CFLAGS = -std=c99 -O3 -D_BSD_SOURCE -D_XOPEN_SOURCE=600

LIBS_ROOT = $(PWD)/../../..
SRC_DIR = $(PWD)
COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(COMMONS_DIR)
LIBS = -lcheck -lcprops -lcurl -Wl,-Bsymbolic-functions


all: string_utils.o region.o

region.o: region.c
	$(CC) $(CFLAGS) -c region.c $(INCLUDES) $(LIBS)
	
string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

clean:
	rm region.o

