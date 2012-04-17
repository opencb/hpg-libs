CC = gcc
CFLAGS = -std=c99 -O3

ROOT = ~/appl/bioinfo-c

SRC_DIR = ..
CONTAINERS_DIR = $(ROOT)/containers
COMMONS_DIR = $(ROOT)/commons
REGION_DIR = $(ROOT)/bioinfo-data/features/region

INCLUDES = -I $(SRC_DIR) -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lpthread -lcheck -lcprops -fopenmp


SRC_OBJECTS = $(SRC_DIR)/vcf_filters.c $(SRC_DIR)/vcf_file.c $(SRC_DIR)/vcf_read.c $(SRC_DIR)/vcf_reader.c $(SRC_DIR)/vcf_write.c $(SRC_DIR)/vcf_batch.c $(SRC_DIR)/util.c
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(REGION_DIR)/region.o
MISC_OBJECTS = list.o string_utils.o

FILTER_OBJECTS = test_vcf_filters.c $(SRC_OBJECTS) $(REGION_TABLE_OBJECTS) $(MISC_OBJECTS)



all: list.o string_utils.o region.o region_table.o test_filter

test_filter: $(FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)

list.o:
	cd $(CONTAINERS_DIR) && make list.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && make region_table.o

clean:
	rm list.o
	rm string_utils.o
	rm region.o
	rm region_table.o
	rm test_filter

