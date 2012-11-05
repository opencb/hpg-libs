
ifeq ($(compiler), intel)
	CC = /opt/intel/bin/icc 
else
	CC = gcc
	compiler=gcc
endif

# -xSSE4.2 -msse4.2 -march=native 

#CC = gcc

CFLAGS = -std=c99 -O3 -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION -g
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION

# Main folders
SRC_DIR = $(PWD)/src
INC_DIR = $(PWD)/include
LIB_DIR = $(PWD)/lib
BIN_DIR = $(PWD)/bin

BIOINFO_LIBS_DIR = $(LIB_DIR)/bioinfo-libs
COMMON_LIBS_DIR = $(LIB_DIR)/common-libs

CONTAINERS_DIR = $(COMMON_LIBS_DIR)/containers
COMMONS_DIR = $(COMMON_LIBS_DIR)/commons
BIOFORMATS_DIR = $(BIOINFO_LIBS_DIR)/bioformats
ALIGNERS_DIR = $(BIOINFO_LIBS_DIR)/aligners

# Include and lib folders
INCLUDES = -I . -I ./include -I $(LIB_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR) -I $(INC_DIR) -I /usr/include/libxml2 -I/usr/local/include
LIBS = -L$(LIB_DIR) -L/usr/lib/x86_64-linux-gnu -Wl,-Bsymbolic-functions -lcprops -fopenmp -largtable2 -lconfig -lbam -lcurl -lm -lz -lxml2

#CUDA_HOME = /usr/local/cuda
CUDA_HOME = /opt/cuda/4.1/cuda
CUDA_LIBS = -L $(CUDA_HOME)/lib64 -lcudart
#CUDA_LIBS = -L /usr/local/cuda/lib64 -lcudart

# -largtable2 -lconfig 

# Object file dependencies
MISC_OBJS = $(COMMONS_DIR)/*.o $(CONTAINERS_DIR)/*.o $(BIOFORMATS_DIR)/fastq/*.o $(BIOFORMATS_DIR)/bam-sam/*.o $(ALIGNERS_DIR)/bwt/*.o $(ALIGNERS_DIR)/sw/*.o 

# Project source files
HPG_ALIGNER_FILES = sw.c rna_server.c rna_splice.c bwt_server_cpu.c region_seeker.c cal_seeker.c pair_server.c sw_server.c batch_writer.c batch_aligner.c buffers.c timing.c genome.c statistics.c options.c

# Project object files
HPG_ALIGNER_OBJS = sw.o rna_server.o rna_splice.o bwt_server_cpu.o region_seeker.o cal_seeker.o pair_server.o sw_server.o batch_writer.o batch_aligner.o buffers.o timing.o genome.o statistics.o options.c

ALL_OBJS = $(HPG_ALIGNER_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies hpg-aligner

hpg-aligner-gpu: compile-dependencies compile-dependencies-gpu  
	cd $(SRC_DIR) &&                                                         \
	$(CC) $(CFLAGS) -DHPG_GPU -c main.c $(HPG_ALIGNER_FILES) $(INCLUDES) $(LIBS) $(CUDA_LIBS) &&    \
	$(CC) $(CFLAGS) -DHPG_GPU -o $(BIN_DIR)/$@ main.o $(ALL_OBJS) $(INCLUDES) $(LIBS) $(CUDA_LIBS)

hpg-aligner: compile-dependencies
	cd $(SRC_DIR) &&                                                         \
	$(CC) $(CFLAGS) -c main.c $(HPG_ALIGNER_FILES) $(INCLUDES) $(LIBS) &&    \
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@ main.o $(ALL_OBJS) $(INCLUDES) $(LIBS)

preprocess-dna: 
	cd $(SRC_DIR) &&                                                         \
	$(CC) $(CFLAGS) -c preprocess_dna.c genome.c $(INCLUDES) $(LIBS) &&    \
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@  preprocess_dna.o genome.o $(INCLUDES) $(LIBS)

compile-dependencies: bam-dependencies
	cd $(COMMONS_DIR) && make compiler=$(compiler) &&           \
	cd $(CONTAINERS_DIR) && make compiler=$(compiler) &&        \
	cd $(ALIGNERS_DIR)/bwt && make compiler=$(compiler) &&      \
	cd $(ALIGNERS_DIR)/sw && make compiler=$(compiler) &&       \
	cd $(BIOFORMATS_DIR)/fastq && make compiler=$(compiler)

compile-dependencies-gpu:
	cd $(ALIGNERS_DIR)/bwt && make test-search-gpu

bam-dependencies:
	cd $(BIOFORMATS_DIR)/bam-sam &&  \
        $(CC) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/bam-sam/alignment.o $(BIOFORMATS_DIR)/bam-sam/alignment.c $(INCLUDES) $(LIBS) && \
        $(CC) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/bam-sam/bam_file.o $(BIOFORMATS_DIR)/bam-sam/bam_file.c $(INCLUDES) $(LIBS)

clean:
	-rm -f $(SRC_DIR)/*~ $(SRC_DIR)/\#*\# $(SRC_DIR)/*.o 
	-rm -f $(COMMONS_DIR)/*.o
	-rm -f $(CONTAINERS_DIR)/*.o
	-rm -f $(BIOFORMATS_DIR)/fastq/*.o
	-rm -f $(BIOFORMATS_DIR)/bam-sam/alignment.o
	-rm -f $(ALIGNERS_DIR)/bwt/*.o
	-rm -f $(ALIGNERS_DIR)/sw/*.o
	-rm -f $(BIN_DIR)/hpg-aligner-gpu
	-rm -f $(BIN_DIR)/hpg-aligner
























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
