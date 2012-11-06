SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp $(SIMD) -D_GNU_SOURCE

SW_DIR = aligners/sw
BWT_DIR = aligners/bwt
COMMON_DIR = ~/appl/bioinfo-c/libs/common-libs/commons
CONTAINERS_LIBS_DIR = ~/appl/bioinfo-c/libs/common-libs/containers
BIOFORMATS_LIBS_DIR = ~/appl/bioinfo-c/libs/bioinfo-libs/bioformats

FASTQ_DIR = $(BIOFORMATS_LIBS_DIR)/fastq
BAM_DIR  = $(BIOFORMATS_LIBS_DIR)/bam-sam

MISC_OBJECTS = $(FASTQ_DIR)/fastq_file.o $(FASTQ_DIR)/fastq_read.o $(FASTQ_DIR)/fastq_batch.o $(FASTQ_DIR)/fastq_batch_reader.o $(BAM_DIR)/alignment.o $(CONTAINERS_LIBS_DIR)/array_list.o $(CONTAINERS_LIBS_DIR)/list.o $(COMMON_DIR)/log.o $(COMMON_DIR)/system_utils.o $(COMMON_DIR)/string_utils.o


INCLUDES = -I . -I $(COMMON_DIR)
LIBS = -L $(SAMTOOLS_DIR) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm -lbam -lz

#SW_FILES = $(SW_DIR)/macros.c $(SW_DIR)/sse.c $(SW_DIR)/emboss.c $(SW_DIR)/smith_waterman.c 
SW_FILES = $(SW_DIR)/*.c 
BWT_FILES = $(BWT_DIR)/*.c $(COMMON_DIR)/string_utils.c

BIOFORMATS_FILES = $(BIOFORMATS_LIBS_DIR)/bam-sam/*.c $(BIOFORMATS_LIBS_DIR)/family/*.c $(BIOFORMATS_LIBS_DIR)/fastq/*.c $(BIOFORMATS_LIBS_DIR)/features/*.c 

MAIN_OBJECTS = *.o 

#$(BWT_DIR)/*.c $(BWT_DIR)/*.o

all: build_aligners

build_aligners:
	$(CC) $(CFLAGS) -fPIC -c $(SW_FILES) $(BWT_FILES) $(INCLUDES)
	$(CC) -shared -o libaligners.so $(MAIN_OBJECTS)
	make clean

compile-dependencies:
	cd $(CONTAINERS_LIBS_DIR) && make compiler=$(compiler) compile && \
	cd $(COMMON_LIBS_DIR) && make compiler=$(compiler) compile && \
	cd $(FASTQ_DIR) && make compiler=$(compiler) compile && \
	cd $(BAM_DIR) && make compiler=$(compiler) compile

clean: 
	rm *.o

## Aligners target
build_sw:
	@echo "Invoking Smith-Waterman makefile to build"
	(cd $(SW_DIR); $(MAKE))

## bioformats target 
build_bwt:
	@echo "Invoking Burrows-Wheeler-transform to build"



