CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Directories' path
LIBS_ROOT = ~/appl/bioinfo-c/libs
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
REGION_DIR = $(LIBS_ROOT)/bioformats/features/region

# Inclusion paths and libraries
INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(LIBS_ROOT)/bioformats/vcf/
LIBS = -lcprops -fopenmp

# Source and object files
MAIN_OBJECTS = ped_file.c ped_read.c ped_reader.c ped_write.c ped_batch.c
MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(REGION_DIR)/region.o
PED_BATCH_OBJECTS = main_batch.c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS)


# Targets
all: compile_dependencies ragel compile main_batch

compile:
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

ragel:
	ragel -o ped_reader.c -C -G2 ped.ragel

ragel_graphviz:
	ragel -o ped_ragel.dot -V -C -G2 ped.ragel

main_batch: compile_dependencies compile ragel $(PED_BATCH_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(PED_BATCH_OBJECTS) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(PED_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

compile_dependencies:
	cd $(COMMONS_DIR) && make file_utils.o && \
	cd $(CONTAINERS_DIR) && make list.o  && \
	cd $(COMMONS_DIR) && make log.o && \
	cd $(REGION_DIR) && make region.o && \
	cd $(CONTAINERS_DIR) && make region_table.o && \
	cd $(COMMONS_DIR) && make string_utils.o	

clean:
	rm -f *.o
	rm main_batch

