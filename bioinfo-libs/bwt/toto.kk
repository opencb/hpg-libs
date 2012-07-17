In file included from ../../../libs/containers/array_list.h:9,
                 from bwt.h:9,
                 from bwt.c:1:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from bwt.h:12,
                 from bwt.c:1:
../../../samtools-0.1.18/bgzf.h: In function 'bgzf_getc':
../../../samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function 'ftello'
In file included from bwt.h:13,
                 from bwt.c:1:
BW_io.h:162:1: warning: "tic" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from bwt.h:12,
                 from bwt.c:1:
../../../libs/commons/commons.h:67:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from bwt.c:1:
BW_io.h:167:1: warning: "toc" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from bwt.h:12,
                 from bwt.c:1:
../../../libs/commons/commons.h:72:1: warning: this is the location of the previous definition
bwt.c: In function 'read_cals_new':
bwt.c:150: warning: assignment discards qualifiers from pointer target type
bwt.c: In function 'read_cals_free':
bwt.c:157: warning: passing argument 2 of 'array_list_free' from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected 'void (*)(void *)' but argument is of type 'void (*)(struct cal_t *)'
bwt.c: In function 'bwt_index_new':
bwt.c:212: warning: implicit declaration of function 'strdup'
bwt.c:212: warning: assignment makes pointer from integer without a cast
bwt.c: In function 'bwt_map_exact_batch':
bwt.c:607: warning: passing argument 2 of 'array_list_free' from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected 'void (*)(void *)' but argument is of type 'void (*)(struct fastq_read_t *)'
bwt.c: In function 'bwt_map_inexact_batch':
bwt.c:1273: warning: passing argument 2 of 'array_list_free' from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected 'void (*)(void *)' but argument is of type 'void (*)(struct fastq_read_t *)'
bwt.c: In function 'bwt_find_cals_from_batch':
bwt.c:1697: warning: passing argument 1 of 'bwt_generate_cal_list' from incompatible pointer type
bwt.c:1403: note: expected 'struct array_list_t *' but argument is of type 'struct array_list_t **'
bwt.c:1728: warning: passing argument 1 of 'array_list_free' from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected 'struct array_list_t *' but argument is of type 'struct array_list_t **'
BW_io.c: In function 'load_exome_file':
BW_io.c:651: warning: initialization makes integer from pointer without a cast
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/file_utils.h:13,
                 from ../../../libs/bioformats/fastq/fastq_file.h:8,
                 from ../../../libs/bioformats/fastq/fastq_file.c:7:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../libs/bioformats/fastq/fastq_batch.c:7:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.h:11,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.c:11:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
../../../libs/bioformats/fastq/fastq_batch_reader.c: In function 'fastq_batch_reader':
../../../libs/bioformats/fastq/fastq_batch_reader.c:40: warning: implicit declaration of function 'omp_get_thread_num'
In file included from ../../../samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
../../../samtools-0.1.18/bgzf.h: In function 'bgzf_getc':
../../../samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function 'ftello'
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../libs/bioformats/bam-sam/alignment.c:3:
../../../libs/bioformats/bam-sam/bam_file.h:7:14: warning: extra tokens at end of #ifdef directive
../../../libs/bioformats/bam-sam/bam_file.h:211:14: warning: extra tokens at end of #ifdef directive
../../../libs/bioformats/bam-sam/alignment.c: In function 'convert_to_quality_string_length':
../../../libs/bioformats/bam-sam/alignment.c:417: warning: function returns address of local variable
../../../libs/bioformats/bam-sam/alignment.c: In function 'convert_to_quality_string':
../../../libs/bioformats/bam-sam/alignment.c:428: warning: function returns address of local variable
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/log.c:1:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/commons/system_utils.h:10,
                 from ../../../libs/commons/system_utils.c:4:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
../../../libs/commons/system_utils.c: In function 'get_optimal_cpu_num_threads':
../../../libs/commons/system_utils.c:69: warning: implicit declaration of function 'popen'
../../../libs/commons/system_utils.c:69: warning: assignment makes pointer from integer without a cast
../../../libs/commons/system_utils.c:82: warning: implicit declaration of function 'pclose'
In file included from ../../../libs/commons/string_utils.c:1:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
../../../libs/commons/string_utils.c: In function 'splitn':
../../../libs/commons/string_utils.c:342: warning: implicit declaration of function 'strsep'
../../../libs/commons/string_utils.c:342: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c:343: warning: implicit declaration of function 'strdup'
../../../libs/commons/string_utils.c:343: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c: In function 'get_to_first_blank':
../../../libs/commons/string_utils.c:441: warning: return makes integer from pointer without a cast
In file included from ../../../libs/containers/array_list.h:9,
                 from ../../../libs/containers/array_list.c:4:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
../../../libs/containers/array_list.c: In function 'array_list_swap':
../../../libs/containers/array_list.c:383: warning: initialization makes integer from pointer without a cast
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/containers/list.c:5:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../libs/containers/array_list.h:9,
                 from test1.h:9,
                 from test1.c:1:
../../../libs/commons/string_utils.h:18:1: warning: "CIGAR_PADDING" redefined
../../../libs/commons/string_utils.h:16:1: warning: this is the location of the previous definition
In file included from ../../../samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from test1.h:12,
                 from test1.c:1:
../../../samtools-0.1.18/bgzf.h: In function 'bgzf_getc':
../../../samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function 'ftello'
In file included from bwt.h:13,
                 from test1.h:13,
                 from test1.c:1:
BW_io.h:162:1: warning: "tic" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from test1.h:12,
                 from test1.c:1:
../../../libs/commons/commons.h:67:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from test1.h:13,
                 from test1.c:1:
BW_io.h:167:1: warning: "toc" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from test1.h:12,
                 from test1.c:1:
../../../libs/commons/commons.h:72:1: warning: this is the location of the previous definition
test1.c: In function 'main':
test1.c:62: warning: implicit declaration of function 'bwt_map_exact_seq'
test1.c:65: warning: implicit declaration of function 'bwt_map_inexact_seq'
