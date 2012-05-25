CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Directories' paths
# LIBS_ROOT = ~/appl/bioinfo-c/libs
LIBS_ROOT = $(PWD)/../..
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
BIOFORMATS_DIR = $(LIBS_ROOT)/bioformats
REGION_DIR = $(BIOFORMATS_DIR)/features/region

# Include path and libraries
INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(BIOFORMATS_DIR)/gff
LIBS = -lcprops -fopenmp

# Source and object files
MAIN_OBJECTS = vcf_file.c vcf_read.c vcf_reader.c vcf_write.c vcf_batch.c vcf_filters.c vcf_stats.c vcf_util.c
MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o $(REGION_DIR)/region.o
GFF_OBJECTS = $(BIOFORMATS_DIR)/gff/*.o
VCF_OBJECTS = ./vcf_*.o

VCF_BASIC_OBJECTS = main_basic.c $(MAIN_OBJECTS) $(MISC_OBJECTS)
VCF_BATCH_OBJECTS = $(VCF_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS)
VCF_FILTER_OBJECTS = main_filter.c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS)


# Targets
all: compile_dependencies ragel compile main_batch main_filter

compile:
# 	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

ragel:
	ragel -o vcf_reader.c -C -G2 vcf.ragel

ragel_graphviz:
	ragel -o vcf_ragel.dot -V -C -G2 vcf.ragel

main_batch: compile_dependencies compile $(VCF_BATCH_OBJECTS) 
# 	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ main_batch.c $(VCF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ main_batch.c $(VCF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

main_filter: compile_dependencies ragel compile $(VCF_FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(VCF_FILTER_OBJECTS) $(INCLUDES) $(LIBS)
# 	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(VCF_FILTER_OBJECTS) $(INCLUDES) $(LIBS)

main_basic: $(VCF_BASIC_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(VCF_BASIC_OBJECTS) $(INCLUDES) $(LIBS)

main_sizeof: main_sizeof_structs.c
	$(CC) $(CFLAGS) -o $@ main_sizeof_structs.c $(INCLUDES)

compile_dependencies:
	cd $(CONTAINERS_DIR) && make list.o  && \
	cd $(COMMONS_DIR) && make log.o && \
	cd $(REGION_DIR) && make region.o && \
	cd $(CONTAINERS_DIR) && make region_table.o && \
	cd $(CONTAINERS_DIR) && make region_table_utils.o && \
	cd $(COMMONS_DIR) && make string_utils.o && \
	cd $(COMMONS_DIR) && make file_utils.o && \
	cd $(BIOFORMATS_DIR)/gff && make compile

clean:
	rm -f *.o
	rm main_batch
	rm main_filter

