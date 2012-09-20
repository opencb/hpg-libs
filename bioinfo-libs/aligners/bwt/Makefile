
LIBS_DIR = ../../../libs
COMMONS_DIR = $(LIBS_DIR)/commons
COMMONS_CUDA_DIR = $(LIBS_DIR)/commons-cuda
CONTAINERS_DIR = $(LIBS_DIR)/containers
FASTQ_DIR = $(LIBS_DIR)/bioformats/fastq
SAM_DIR = $(LIBS_DIR)/bioformats/bam-sam
#SAMTOOLS_DIR = ../../../hpg-aligner/samtools-0.1.18/
#SAMTOOLS_DIR = /home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18
SAMTOOLS_DIR = ../../../samtools-0.1.18/
CPROPS_DIR = ../../../libcprops-0.1.12/

LIB = -L$(SAM_DIR) -L$(SAMTOOLS_DIR) -L$(CPROPS_DIR) -lpthread -lcprops -lbam -lm -lz -fopenmp 

ALL = test-preprocess test-search-cpu

CC = gcc
CFLAGS =  -O3 -fopenmp -std=c99 -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION -DVERBOSE
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

NVCC = nvcc
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler " -fopenmp"
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12 -DVECTOR_O_COMPRESSION
#NVCCFLAGS = -g -G -Xptxas -v

CINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(CPROPS_DIR) -I$(SAMTOOLS_DIR) -I$(SAM_DIR) 
CUINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(CPROPS_DIR) -I$(SAMTOOLS_DIR) -I$(SAM_DIR) -I/opt/cuda/include

FASTQ-SRCS = $(FASTQ_DIR)/fastq_file.c $(FASTQ_DIR)/fastq_read.c $(FASTQ_DIR)/fastq_batch.c $(FASTQ_DIR)/fastq_batch_reader.c
SAM-SRCS = $(SAM_DIR)/alignment.c
CONTAINERS-SRCS = $(CONTAINERS_DIR)/array_list.c $(CONTAINERS_DIR)/list.c
COMMONS-SRCS = $(COMMONS_DIR)/log.c $(COMMONS_DIR)/system_utils.c $(COMMONS_DIR)/string_utils.c


SRCS = bwt.c BW_io.c BW_search.c BW_preprocess.c $(FASTQ-SRCS) $(SAM-SRCS) $(COMMONS-SRCS) $(CONTAINERS-SRCS)
OBJS = *.o 

all: $(ALL)

test-preprocess: objects test-preprocess.c
	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test-preprocess.c -o test-preprocess $(LIB) 

test-search-cpu: objects test-search-cpu.c
	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test-search-cpu.c -o test-search-cpu $(LIB) 

#test1-cpu: objects test1.c test1.h
#	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test1.c -o test1-cpu $(LIB) 

#test2-cpu: objects test2.c test2.h
#	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test2.c -o test2-cpu $(LIB) 

#test3-cpu: objects test3.c test3.h
#	$(CC) $(CFLAGS) $(CINCLUDES) $(OBJS) test3.c -o test3-cpu $(LIB) 

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

