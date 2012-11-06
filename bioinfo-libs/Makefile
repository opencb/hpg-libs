SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp $(SIMD) -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION

SW_DIR = aligners/sw
BWT_DIR = aligners/bwt

COMMON_LIB_DIR = ~/appl/bioinfo-c/libs/common-libs/
COMMON_DIR = $(COMMON_LIB_DIR)/commons
CONTAINERS_DIR = $(COMMON_LIB_DIR)/containers
BIOFORMATS_LIBS_DIR = ~/appl/bioinfo-c/libs/bioinfo-libs/bioformats
SAMTOOLS_DIR = ~/appl/bioinfo-c/libs/ext/samtools-0.1.18
BWT_DIR = ~/appl/bioinfo-c/libs/bioinfo-libs/aligners/bwt

FASTQ_DIR = $(BIOFORMATS_LIBS_DIR)/fastq
BAM_DIR  = $(BIOFORMATS_LIBS_DIR)/bam-sam

MISC_OBJECTS = $(FASTQ_DIR)/fastq_file.o $(FASTQ_DIR)/fastq_read.o $(FASTQ_DIR)/fastq_batch.o $(FASTQ_DIR)/fastq_batch_reader.o $(BAM_DIR)/alignment.o $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMON_DIR)/log.o $(COMMON_DIR)/system_utils.o $(COMMON_DIR)/string_utils.o


INCLUDES = -I . -I $(COMMON_LIB_DIR) -I $(SAMTOOLS_DIR) -I $(BWT_DIR)
LIBS = -L $(SAMTOOLS_DIR) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm -lbam -lz

ALIGNER_FILES = $(SW_DIR)/*.c $(BWT_DIR)/*.c  
BIOFORMATS_FILES = $(BIOFORMATS_LIBS_DIR)/bam-sam/alignment.c $(BIOFORMATS_LIBS_DIR)/family/*.c $(BIOFORMATS_LIBS_DIR)/fastq/*.c $(BIOFORMATS_LIBS_DIR)/features/region/*.c $(BIOFORMATS_LIBS_DIR)/vcf/*.c $(BIOFORMATS_LIBS_DIR)/gff/*.c $(BIOFORMATS_LIBS_DIR)/ped/*.c

MAIN_OBJECTS = *.o

#$(BWT_DIR)/*.c $(BWT_DIR)/*.o

all: build_aligners

build_aligners: compile-dependencies
	$(CC) $(CFLAGS) -fPIC -c $(ALIGNER_FILES) $(BIOFORMATS_FILES) $(INCLUDES)
	$(CC) -shared -o libbioinfo.so $(MAIN_OBJECTS)
	make clean

compile-dependencies:
	cd $(CONTAINERS_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(COMMON_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(FASTQ_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(BAM_DIR) && make COMPILER=$(COMPILER) compile

clean: 
	rm *.o
