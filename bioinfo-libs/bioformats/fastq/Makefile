CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE -g
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT)/common-libs
LIBS = -fopenmp -Wl,-Bsymbolic-functions

# Source and object files
MAIN_OBJECTS = fastq_batch.c fastq_batch_reader.c fastq_file.c fastq_read.c
#MISC_OBJECTS = $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o

# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make && \
	cd $(COMMONS_DIR) && make

clean:
	rm -f *.o
