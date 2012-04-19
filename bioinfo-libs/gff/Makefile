CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = ../..
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
REGION_DIR = $(LIBS_ROOT)/bioformats/features/region

INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lcprops -fopenmp

MAIN_OBJECTS = vcf_file.c vcf_read.c vcf_reader.c vcf_write.c vcf_batch.c $(CONTAINERS_DIR)/list.o
VCF_BASIC_OBJECTS = main_basic.c $(MAIN_OBJECTS)
VCF_BATCH_OBJECTS = main_batch.c $(MAIN_OBJECTS)

REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(REGION_DIR)/region.o $(COMMONS_DIR)/string_utils.o
VCF_FILTER_OBJECTS = main_filter.c $(MAIN_OBJECTS) vcf_filters.c util.c $(REGION_TABLE_OBJECTS)



all: list.o string_utils.o region.o region_table.o ragel main_batch main_filter

ragel:
	ragel -o vcf_reader.c -C -G2 vcf.ragel

ragel_graphviz:
	ragel -o vcf_ragel.dot -V -C -G2 vcf.ragel

main_batch: $(VCF_BATCH_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(VCF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

main_filter: $(VCF_FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(VCF_FILTER_OBJECTS) $(INCLUDES) $(LIBS)
# 	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(VCF_FILTER_OBJECTS) $(INCLUDES) $(LIBS)

main_basic: $(VCF_BASIC_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(VCF_BASIC_OBJECTS) $(INCLUDES) $(LIBS)

main_sizeof: main_sizeof_structs.c
	$(CC) $(CFLAGS) -o $@ main_sizeof_structs.c $(INCLUDES)

list.o:
	cd $(CONTAINERS_DIR) && make list.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && make region_table.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o


clean:
	rm main_batch
	rm main_filter

