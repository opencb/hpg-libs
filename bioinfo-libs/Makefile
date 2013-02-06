SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp $(SIMD) -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION

## This param must point to samtools project, 'make' is assumed
SAMTOOLS_DIR = ~/appl/bioinfo-c/libs/ext/samtools-0.1.18

SW_DIR = aligners/sw
BWT_DIR = aligners/bwt

COMMON_LIB_DIR = ~/appl/bioinfo-c/libs/common-libs/
COMMON_DIR = $(COMMON_LIB_DIR)/commons
CONTAINERS_DIR = $(COMMON_LIB_DIR)/containers
BIOFORMATS_LIBS_DIR = ~/appl/bioinfo-c/libs/bioinfo-libs/bioformats

FASTQ_DIR = $(BIOFORMATS_LIBS_DIR)/fastq
BAM_DIR  = $(BIOFORMATS_LIBS_DIR)/bam-sam
#MISC_OBJECTS = $(FASTQ_DIR)/fastq_file.o $(FASTQ_DIR)/fastq_read.o $(FASTQ_DIR)/fastq_batch.o $(FASTQ_DIR)/fastq_batch_reader.o $(BAM_DIR)/alignment.o $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMON_DIR)/log.o $(COMMON_DIR)/system_utils.o $(COMMON_DIR)/string_utils.o

INCLUDES = -I . -I $(COMMON_LIB_DIR) -I $(SAMTOOLS_DIR) -I $(BWT_DIR)
LIBS = -L $(SAMTOOLS_DIR) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm -lbam -lz

ALIGNER_FILES = $(SW_DIR)/*.c $(BWT_DIR)/*.c  
BIOFORMATS_FILES = $(BIOFORMATS_LIBS_DIR)/bam-sam/alignment.c $(BIOFORMATS_LIBS_DIR)/family/*.c $(BIOFORMATS_LIBS_DIR)/fastq/*.c $(BIOFORMATS_LIBS_DIR)/features/region/*.c $(BIOFORMATS_LIBS_DIR)/vcf/*.c $(BIOFORMATS_LIBS_DIR)/gff/*.c $(BIOFORMATS_LIBS_DIR)/ped/*.c

MAIN_OBJECTS = *.o

all: build_aligners


build_aligners: compile-dependencies
	$(CC) $(CFLAGS) -fPIC -c $(ALIGNER_FILES) $(BIOFORMATS_FILES) $(INCLUDES)
	$(CC) -shared -o libbioinfo.so $(MAIN_OBJECTS)
	rm *.o

compile-dependencies:
	cd $(CONTAINERS_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(COMMON_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(FASTQ_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(BAM_DIR) && make COMPILER=$(COMPILER) compile

clean: 
	rm *.o *.so









#BIN = ../bin
#LIB = -L./samtools-0.1.18/ -L/usr/local/cuda/lib64/ -L./libcprops-0.1.12/ -lpthread -lcprops -lbam -lz -lm -L./libs/bioformats/bam-sam/

#ALL = main-cpu

#CC = gcc
#CFLAGS =  -g -O3 -Wall -fopenmp -std=c99 -D_GNU_SOURCE -DVECTOR_O_COMPRESSION -DVERBOSE 

#CINCLUDES = -I. -Ilibs/aligners/bwt/ -Ilibs/aligners/sw/ -I./samtools-0.1.18/ -I. -Ilibs/commons-cuda/ -Ilibs/commons/ -Ilibs/containers/ -Ilibs/bioformats/bam-sam/ -Ilibs/bioformats/fastq/ -I/opt/cuda/include/ -I./samtools-0.1.18/ -I./libcprops-0.1.12/

#CUINCLUDES = -I. -Ilibs/aligners/bwt/ -I./samtools-0.1.18/ -I. -Ilibs/commons-cuda/ -Ilibs/commons/ -Ilibs/containers/ -Ilibs/bioformats/bam-sam/ -Ilibs/bioformats/fastq/ -I/opt/cuda/include/ -I./samtools-0.1.18/ -I./libcprops-0.1.12/

#NVCC = nvcc
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12


#COMMONS-ALIGNER-SOURCE = libs/bioformats/fastq/fastq_batch.c libs/bioformats/bam-sam/alignment.c ./main.c libs/bioformats/fastq/fastq_batch_reader.c ./buffers.c ./timing.c ./genome.c libs/aligners/bwt/bwt.c libs/aligners/sw/macros.c libs/aligners/sw/sse.c libs/aligners/sw/smith_waterman.c libs/containers/array_list.c libs/commons/log.c libs/containers/list.c  ./bwt_server_cpu.c libs/aligners/bwt/BW_io.c libs/aligners/bwt/BW_search.c libs/commons/string_utils.c libs/bioformats/fastq/fastq_file.c libs/bioformats/fastq/fastq_read.c ./region_seeker.c ./cal_seeker.c ./sw_server.c ./pair_server.c ./batch_writer.c ./batch_aligner.c libs/bioformats/bam-sam/bam_file.c ./statistics.c


#CPU-ALIGNER-OBJ = ./fastq_batch.o ./alignment.o ./main.o ./fastq_batch_reader.o ./buffers.o ./timing.o ./genome.o ./bwt.o ./macros.o ./sse.o ./smith_waterman.o ./array_list.o ./log.o ./list.o ./bwt_server_cpu.o ./BW_io.o ./BW_search.o ./string_utils.o ./fastq_file.o ./fastq_read.o ./region_seeker.o ./cal_seeker.o ./sw_server.o ./pair_server.o ./batch_writer.o ./batch_aligner.o ./bam_file.o ./statistics.o

#GPU-ALIGNER-OBJ = ./bwt_server_cuda.o

#all: $(ALL)

#sequence: genome.c sequence.c
#	$(CC) $(CFLAGS) genome.c sequence.c -o sequence

#main-cuda: cpu-objects cuda-objects
#	$(CC) $(CFLAGS) -DCUDA_VERSION $(CPU-ALIGNER-OBJ) $(GPU-ALIGNER-OBJ) -o main $(LIB)  -I./libcprops-0.1.12 -lcudart

#main-cpu: cpu-objects
#	$(CC) $(CFLAGS) $(CPU-ALIGNER-OBJ) -o main-cpu $(LIB)

#cpu-objects: $(COMMONS-ALIGNER-SOURCE)
#	$(CC) $(CFLAGS) $(CINCLUDES) -c $^

#cuda-objects: ./bwt_server_cuda.cu
#	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION $(CUINCLUDES) -c $^

##############################################################################3
#clean-bam-gpu-tools:
#	rm -f *~ \#*\# $(LIB)/*.o $(BIN)/bam-gpu-tools

#clean:
#	rm -f *~ \#*\# *.o $(BIN)/* $(ALL)
