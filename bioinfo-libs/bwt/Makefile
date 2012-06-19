LIB = -L/usr/local/cuda/lib64/

LIBS_DIR = ../../../libs
COMMONS_DIR = $(LIBS_DIR)/commons
COMMONS_CUDA_DIR = $(LIBS_DIR)/commons-cuda
CONTAINERS_DIR = $(LIBS_DIR)/containers
FASTQ_DIR = $(LIBS_DIR)/bioformats/fastq
SAM_DIR = $(LIBS_DIR)/bioformats/bam-sam
SAMTOOLS_DIR = /home/jtarraga/softs/samtools-0.1.17

ALL = test1-cpu

CC = gcc
CFLAGS =  -O3 -fopenmp -std=c99 -DVECTOR_O_COMPRESSION -DVERBOSE
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

NVCC = nvcc
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler " -fopenmp"
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12 -DVECTOR_O_COMPRESSION
#NVCCFLAGS = -g -G -Xptxas -v

CINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(SAM_DIR) -I$(SAMTOOLS_DIR) 
CUINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(SAM_DIR) I$(SAMTOOLS_DIR) -I/opt/cuda/include

FASTQ-SRCS = $(FASTQ_DIR)/fastq_file.c $(FASTQ_DIR)/fastq_read.c $(FASTQ_DIR)/fastq_batch.c

COMMONS-SRCS = $(COMMONS_DIR)/log.c $(COMMONS_DIR)/system_utils.c $(COMMONS_DIR)/string_utils.c
#CONTAINERS-SRCS = $(CONTAINERS_DIR)/list.c

#ALIGNER-SRCS = main.c alignment_batch.c fastq_batch_reader.c result_server.c timing.c BW_io.c BW_search.c
#GPU-ALIGNER-OBJ = *.o #main.o fastq_batch_reader.o timing.c $(FASTQ_DIR)/*.o $(CONTAINERS_DIR)/list.o 

SRCS = burrows_wheeler.c BW_io.c BW_search.c $(FASTQ-SRCS) $(COMMONS-SRCS)
OBJS = *.o #burrows_wheeler.o BW_io.o 

all: $(ALL)

test1-cpu: objects test1.c test1.h
	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test1.c -o test1-cpu $(LIB)

objects: $(SRCS)
	$(CC) $(CFLAGS) $(CINCLUDES) -c $^

#main-cuda: cpu-objects cuda-objects
#        $(CC) $(CFLAGS) -DCUDA_VERSION $(GPU-ALIGNER-OBJ) -o main $(LIB) -lcudart

#main-omp: omp-objects
#        $(CC) $(CFLAGS) $(GPU-ALIGNER-OBJ) -o main $(LIB) -lpthread 

#omp-objects: $(ALIGNER-SRCS) bwt_server_omp.c
#        $(CC) $(CFLAGS) $(CINCLUDES) -c $^

#cpu-objects: $(COMMONS-SRCS) $(CONTAINERS-SRCS) $(FASTQ-SRCS) $(ALIGNER-SRCS) 
#        $(CC) $(CFLAGS) -DCUDA_VERSION $(CINCLUDES) -c $^

#cuda-objects: BW_gpu.cu bwt_context_cuda.cu gpu_context_cuda.cu exact_seeker_cuda.cu #./bwt_server_cuda.cu
#        $(NVCC) $(NVCCFLAGS) -DCUDA_VERSION $(CUINCLUDES) -c $^

#containers-objects:
#       (cd $(CONTAINERS_DIR); $(MAKE))

#sequence: genome.c sequence.c
#        $(CC) $(CFLAGS) genome.c sequence.c -o sequence

###################################################################

clean:
	-rm -f *~ \#*\# *.o $(ALL) main

