BIN = ../bin
LIB = .

ALL = io

CC = gcc
CFLAGS = -Wall -O3 -std=c99
CFLAGS_DEBUG = -Wall -pg

LIB_ROOT = ..
CUDA_COMMONS_DIR = $(LIB_ROOT)/commons-cuda

LIBS = -lcurl -Wl,-Bsymbolic-functions


all: system_utils.o file_utils.o http_utils.o log.o string_utils.o result.o
test: system_utils.o file_utils.o http_utils.o log.o string_utils.o result.o test-utils

system_utils.o: system_utils.h system_utils.c
	$(CC) $(CFLAGS) -DCUDA_VERSION -c system_utils.c -I$(CUDA_COMMONS_DIR)

file_utils.o: file_utils.h file_utils.c string_utils.o
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c file_utils.c

http_utils.o: http_utils.h http_utils.c
	$(CC) $(CFLAGS) -c http_utils.c

log.o: log.h log.c string_utils.o
	$(CC) $(CFLAGS) -DCUDA_VERSION -c log.c

string_utils.o: string_utils.h string_utils.c
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c string_utils.c
#	$(CC) $(CFLAGS_DEBUG) -c string_utils.c

result.o: result.h result.c
	$(CC) $(CFLAGS) -c result.c -I/usr/include/libxml2 -lxml2

test-utils:
	cd test && make

clean:
	-rm -f *.o
