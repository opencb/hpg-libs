#CC = gcc
#CC = /opt/intel/bin/icc
ifeq ($(compiler), intel)
        CC = /opt/intel/bin/icc
else
        CC = gcc
        compiler=gcc
endif

CFLAGS = -std=c99 -O3 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

## This param must point to samtools project, 'make' is assumed
SAMTOOLS_DIR = ~/appl/bioinfo-c/libs/ext/samtools-0.1.18

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I ./include -I $(SAMTOOLS_DIR)
LIBS = -fopenmp -Wl,-Bsymbolic-functions

# Source and object files
MAIN_OBJECTS = alignment.c
#MISC_OBJECTS = $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o

# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make compiler=$(compiler) && \
	cd $(COMMONS_DIR) && make compiler=$(compiler)

clean:
	rm -f *.o
