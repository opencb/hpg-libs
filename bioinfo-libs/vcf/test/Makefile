CC = gcc
CFLAGS = -std=c99 -O3

LIBS_ROOT = ../../..

SRC_DIR = ..
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
REGION_DIR = $(LIBS_ROOT)/bioformats/features/region
GFF_DIR = $(LIBS_ROOT)/bioformats/gff

INCLUDES = -I $(SRC_DIR) -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(GFF_DIR)
LIBS = -lpthread -lcheck -lcprops -fopenmp


SRC_OBJECTS = $(SRC_DIR)/vcf_filters.c $(SRC_DIR)/vcf_file.c $(SRC_DIR)/vcf_read.c $(SRC_DIR)/vcf_reader.c $(SRC_DIR)/vcf_write.c $(SRC_DIR)/vcf_batch.c $(SRC_DIR)/vcf_util.c
GFF_OBJECTS = $(GFF_DIR)/gff_file.c $(GFF_DIR)/gff_read.c $(GFF_DIR)/gff_reader.c $(GFF_DIR)/gff_write.c $(GFF_DIR)/gff_batch.c
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o $(REGION_DIR)/region.o

MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o

# Test executables
FILTER_OBJECTS = test_vcf_filters.c $(SRC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS) $(MISC_OBJECTS)
REGION_PARSE_OBJECTS = test_parse_regions.c $(GFF_OBJECTS) $(REGION_TABLE_OBJECTS) $(MISC_OBJECTS) $(SRC_DIR)/util.c


all: list.o string_utils.o region.o region_table.o test_filter test_parse_region

test_parse_region: $(REGION_PARSE_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(REGION_PARSE_OBJECTS) $(INCLUDES) $(LIBS)

test_filter: $(FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)

list.o:
	cd $(CONTAINERS_DIR) && make list.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && $(MAKE) region_table.o

clean:
	rm list.o
	rm string_utils.o
	rm region.o
	rm region_table.o
	rm test_filter

