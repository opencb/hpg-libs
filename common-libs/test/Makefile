
CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

ROOT = ~/appl/bioinfo-c
SRC_DIR = ..
REGION_DIR = $(ROOT)/bioinfo-data/features/region
COMMONS_DIR = $(ROOT)/libs/commons

INCLUDES = -I $(SRC_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR)
LIBS = -lm -lcprops -lcheck


TEST_REGION_TABLE_OBJS = test_region_table.c $(SRC_DIR)/region_table.o $(COMMONS_DIR)/string_utils.o $(REGION_DIR)/region.o 
TEST_ARRAY_LIST_OBJS = test_array_list.c $(SRC_DIR)/array_list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/log.o

all: region_table.o test_region_table


test_region_table: $(TEST_REGION_TABLE_OBJS)
	$(CC) $(CFLAGS) -g -o $@ $(TEST_REGION_TABLE_OBJS) $(INCLUDES) $(LIBS)

test_array_list: $(TEST_ARRAY_LIST_OBJS)
	$(CC) $(CFLAGS_DEBUG) -g -D_BSD_SOURCE -o $@ $(TEST_ARRAY_LIST_OBJS) $(INCLUDES) $(LIBS)

region_table.o: $(REGION_TABLE)
	cd .. && make region_table.o

array_list.o:
	cd .. && make array_list.o
	
clean:
	test ! -f test_region_table || rm test_region_table
	test ! -f test_array_list || rm test_array_list

