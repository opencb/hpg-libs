CC = gcc
CFLAGS = -std=c99 -O3 -D_BSD_SOURCE -D_XOPEN_SOURCE=600

LIBS_ROOT = $(PWD)/../../..
# SRC_DIR = $(PWD)
# COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I $(PWD) -I $(LIBS_ROOT)
LIBS = -lcprops -lcurl -Wl,-Bsymbolic-functions
LIBS_TEST = -lcheck

MISC_FILES = $(LIBS_ROOT)/commons/string_utils.o

all: region.o

region.o: compile-dependencies region.c
	$(CC) $(CFLAGS) -c region.c $(MISC_FILES) $(INCLUDES) $(LIBS)
	
compile-dependencies:
	cd $(LIBS_ROOT)/commons && make compile

clean:
	rm region.o

