CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
COMMONS_CUDA_DIR = $(LIBS_ROOT)/commons-cuda
HPG_BAM_DIR = $(LIBS_ROOT)/../hpg-bam-tools
HPG_BAM_FROM_SUBMODULE_DIR = $(LIBS_ROOT)/../../hpg-bam-tools
SAMTOOLS_DIR = $(LIBS_ROOT)/bam

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(COMMONS_DIR) -I $(COMMONS_CUDA_DIR) -I $(HPG_BAM_DIR) -I $(HPG_BAM_FROM_SUBMODULE_DIR) -I $(SAMTOOLS_DIR)
LIBS = -fopenmp -Wl,-Bsymbolic-functions

# Source and object files
MAIN_OBJECTS = alignment.c bam_data_batch.c bam_data_batch_list.c bam_file.c bam_reader.c bam_writer.c

# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make && \
	cd $(COMMONS_DIR) && make

clean:
	rm -f *.o
