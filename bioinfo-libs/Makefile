SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp $(SIMD) -D_GNU_SOURCE

## This param must point to samtools project, 'make' is assumed
SAMTOOLS_DIR = ~/appl/bioinfo-c/libs/ext/samtools-0.1.18

SW_DIR = aligners/sw
BWT_DIR = aligners/bwt

## bioinfo-libs dirs
BIOINFO_LIBS_DIR = ~/appl/bioinfo-c/libs/bioinfo-libs
ALIGNERS_LIBS_DIR = $(BIOINFO_LIBS_DIR)/aligners
BIOFORMATS_LIBS_DIR = $(BIOINFO_LIBS_DIR)/bioformats

## common-libs dirs
COMMONS_LIBS_DIR = ~/appl/bioinfo-c/libs/common-libs
COMMONS_DIR = $(COMMONS_LIBS_DIR)/commons
CONTAINERS_LIBS_DIR = $(COMMONS_LIBS_DIR)/containers

FASTQ_DIR = $(BIOFORMATS_LIBS_DIR)/fastq
BAM_DIR  = $(BIOFORMATS_LIBS_DIR)/bam-sam

MISC_OBJECTS = $(FASTQ_DIR)/fastq_file.o $(FASTQ_DIR)/fastq_read.o $(FASTQ_DIR)/fastq_batch.o $(FASTQ_DIR)/fastq_batch_reader.o $(BAM_DIR)/alignment.o $(CONTAINERS_LIBS_DIR)/array_list.o $(CONTAINERS_LIBS_DIR)/list.o $(COMMON_DIR)/log.o $(COMMON_DIR)/system_utils.o $(COMMON_DIR)/string_utils.o


INCLUDES = -I . -I $(BIOINFO_LIBS_DIR) -I $(COMMONS_LIBS_DIR) -I $(SAMTOOLS_DIR)
LIBS = -L $(SAMTOOLS_DIR) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm -lbam -lz

#SW_FILES = $(SW_DIR)/macros.c $(SW_DIR)/sse.c $(SW_DIR)/emboss.c $(SW_DIR)/smith_waterman.c 
SW_FILES = $(SW_DIR)/*.c 
BWT_FILES = $(BWT_DIR)/*.c $(COMMONS_DIR)/string_utils.c

BIOFORMATS_FILES = $(BIOFORMATS_LIBS_DIR)/bam-sam/*.c $(BIOFORMATS_LIBS_DIR)/family/*.c $(BIOFORMATS_LIBS_DIR)/fastq/*.c $(BIOFORMATS_LIBS_DIR)/features/region/*.c 

MAIN_OBJECTS = *.o 

#$(BWT_DIR)/*.c $(BWT_DIR)/*.o

all: build_aligners

build_aligners:
#	$(CC) $(CFLAGS) -fPIC $(INCLUDES) -c $(BIOFORMATS_FILES)
	$(CC) $(CFLAGS) -fPIC -shared $(INCLUDES) $(BIOFORMATS_FILES) -o build/libaligners.so #$(MAIN_OBJECTS)
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



