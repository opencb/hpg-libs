CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

LIBS_ROOT = ../..
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
REGION_DIR = $(LIBS_ROOT)/bioformats/features/region

INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(LIBS_ROOT)/bioformats/vcf/
LIBS = -lcprops -fopenmp

MAIN_OBJECTS = gff_file.c gff_read.c gff_reader.c gff_write.c gff_batch.c
MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(REGION_DIR)/region.o
GFF_BATCH_OBJECTS = main_batch.c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS)


all: log.o file_utils.o list.o string_utils.o region.o region_table.o ragel main_batch

ragel:
	ragel -o gff_reader.c -C -G2 gff.ragel

ragel_graphviz:
	ragel -o gff_ragel.dot -V -C -G2 gff.ragel

main_batch: $(GFF_BATCH_OBJECTS)
#	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(GFF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(GFF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

file_utils.o:
	cd $(COMMONS_DIR) && make file_utils.o

list.o:
	cd $(CONTAINERS_DIR) && make list.o

log.o:
	cd $(COMMONS_DIR) && make log.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && make region_table.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

clean:
	rm main_batch

