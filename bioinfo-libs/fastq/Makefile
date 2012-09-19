CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
LIBS_ROOT = ../..
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT)
LIBS = -lcprops -fopenmp -Wl,-Bsymbolic-functions

INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include
LIBS_STATIC = -L $(LIBS_ROOT) -lcprops -fopenmp

# Source and object files
MAIN_OBJECTS = fastq_batch.c fastq_batch_reader.c fastq_file.c fastq_read.c

# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-static: compile-dependencies-static
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES_STATIC) $(LIBS_STATIC)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(COMMONS_DIR) && make compile && \

compile-dependencies-static:
	cd $(CONTAINERS_DIR) && make compile-static && \
	cd $(COMMONS_DIR) && make compile && \

clean:
	rm -f *.o
