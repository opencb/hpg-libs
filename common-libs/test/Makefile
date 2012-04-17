
CC = gcc
CFLAGS = -std=c99 -O3

ROOT = ~/appl/bioinfo-c
SRC_DIR = ..
REGION_DIR = $(ROOT)/bioinfo-data/features/region
COMMONS_DIR = $(ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lcprops -lcheck


TEST_REGION_TABLE_OBJS = test_region_table.c $(SRC_DIR)/region_table.o $(COMMONS_DIR)/string_utils.o $(REGION_DIR)/region.o 


all: region_table.o test_region_table


test_region_table: $(TEST_REGION_TABLE_OBJS)
	$(CC) $(CFLAGS) -g -o $@ $(TEST_REGION_TABLE_OBJS) $(INCLUDES) $(LIBS)


region_table.o: $(REGION_TABLE)
	cd .. && make region_table.o

clean:
	rm test_region_table

