In file included from ../../../libs/containers/array_list.h:9,
                 from bwt.h:9,
                 from bwt.c:1:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from /home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from bwt.h:12,
                 from bwt.c:1:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function ‘ftello’
In file included from bwt.h:13,
                 from bwt.c:1:
BW_io.h:50:1: warning: "DELETION" redefined
In file included from ../../../libs/containers/array_list.h:9,
                 from bwt.h:9,
                 from bwt.c:1:
../../../libs/commons/string_utils.h:12:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from bwt.c:1:
BW_io.h:52:1: warning: "INSERTION" redefined
In file included from ../../../libs/containers/array_list.h:9,
                 from bwt.h:9,
                 from bwt.c:1:
../../../libs/commons/string_utils.h:11:1: warning: this is the location of the previous definition
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
bwt.c: In function ‘read_cals_new’:
bwt.c:147: warning: assignment discards qualifiers from pointer target type
bwt.c: In function ‘read_cals_free’:
bwt.c:154: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
bwt.c: In function ‘bwt_index_new’:
bwt.c:209: warning: implicit declaration of function ‘strdup’
bwt.c:209: warning: assignment makes pointer from integer without a cast
bwt.c: In function ‘bwt_map_exact_batch’:
bwt.c:591: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct fastq_read_t *)’
bwt.c: In function ‘bwt_map_inexact_seq’:
bwt.c:844: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:850: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:853: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:856: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:863: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:866: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:869: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:876: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:879: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:882: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c: In function ‘bwt_map_inexact_batch’:
bwt.c:1151: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct fastq_read_t *)’
bwt.c: In function ‘bwt_map_inexact_seeds_seq’:
bwt.c:1188: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1188: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1188: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1196: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1196: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1205: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1205: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1205: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1213: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1213: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1218: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c: In function ‘bwt_generate_cal_list’:
bwt.c:1246: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1258: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1258: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1292: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1292: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1296: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1296: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1296: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1296: warning: format ‘%d’ expects type ‘int’, but argument 6 has type ‘size_t’
bwt.c:1301: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1301: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1301: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1301: warning: format ‘%d’ expects type ‘int’, but argument 6 has type ‘size_t’
bwt.c:1305: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1305: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1334: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1334: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1334: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1334: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1344: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1344: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1369: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1369: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1369: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c: In function ‘bwt_find_cals_from_seq’:
bwt.c:1397: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
BW_io.c: In function ‘load_exome_file’:
BW_io.c:651: warning: initialization makes integer from pointer without a cast
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/file_utils.h:13,
                 from ../../../libs/bioformats/fastq/fastq_file.h:8,
                 from ../../../libs/bioformats/fastq/fastq_file.c:7:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
../../../libs/bioformats/fastq/fastq_file.c: In function ‘fastq_fread_num_reads’:
../../../libs/bioformats/fastq/fastq_file.c:53: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:54: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:55: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c: In function ‘fastq_fread_max_size’:
../../../libs/bioformats/fastq/fastq_file.c:89: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:90: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:91: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c: In function ‘fastq_fread_batch_max_size’:
../../../libs/bioformats/fastq/fastq_file.c:130: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:131: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:132: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c: In function ‘fastq_fread_index_positions’:
../../../libs/bioformats/fastq/fastq_file.c:201: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:202: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:203: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from ../../../libs/bioformats/fastq/fastq_batch.c:7:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.h:11,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.c:11:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
../../../libs/bioformats/fastq/fastq_batch_reader.c: In function ‘fastq_batch_reader’:
../../../libs/bioformats/fastq/fastq_batch_reader.c:40: warning: implicit declaration of function ‘omp_get_thread_num’
In file included from /home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function ‘ftello’
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from ../../../libs/bioformats/bam-sam/alignment.c:3:
../../../libs/bioformats/bam-sam/bam_file.h:7:14: warning: extra tokens at end of #ifdef directive
../../../libs/bioformats/bam-sam/bam_file.h:211:14: warning: extra tokens at end of #ifdef directive
../../../libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string_length’:
../../../libs/bioformats/bam-sam/alignment.c:416: warning: function returns address of local variable
../../../libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string’:
../../../libs/bioformats/bam-sam/alignment.c:427: warning: function returns address of local variable
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/log.c:1:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/commons/system_utils.h:10,
                 from ../../../libs/commons/system_utils.c:4:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
../../../libs/commons/system_utils.c: In function ‘get_free_memory’:
../../../libs/commons/system_utils.c:16: warning: ignoring return value of ‘system’, declared with attribute warn_unused_result
../../../libs/commons/system_utils.c:23: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/commons/system_utils.c: In function ‘get_optimal_cpu_num_threads’:
../../../libs/commons/system_utils.c:69: warning: implicit declaration of function ‘popen’
../../../libs/commons/system_utils.c:69: warning: assignment makes pointer from integer without a cast
../../../libs/commons/system_utils.c:82: warning: implicit declaration of function ‘pclose’
../../../libs/commons/system_utils.c:75: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from ../../../libs/commons/string_utils.c:1:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
../../../libs/commons/string_utils.c: In function ‘splitn’:
../../../libs/commons/string_utils.c:342: warning: implicit declaration of function ‘strsep’
../../../libs/commons/string_utils.c:342: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c:343: warning: implicit declaration of function ‘strdup’
../../../libs/commons/string_utils.c:343: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c: In function ‘get_to_first_blank’:
../../../libs/commons/string_utils.c:441: warning: return makes integer from pointer without a cast
In file included from ../../../libs/containers/array_list.h:9,
                 from ../../../libs/containers/array_list.c:4:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/containers/list.c:5:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from ../../../libs/containers/array_list.h:9,
                 from test3.h:9,
                 from test3.c:1:
../../../libs/commons/string_utils.h:17:1: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:1: warning: this is the location of the previous definition
In file included from /home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bam.h:52,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from test3.h:12,
                 from test3.c:1:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
/home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18/bgzf.h:149: warning: implicit declaration of function ‘ftello’
In file included from bwt.h:13,
                 from test3.h:13,
                 from test3.c:1:
BW_io.h:50:1: warning: "DELETION" redefined
In file included from ../../../libs/containers/array_list.h:9,
                 from test3.h:9,
                 from test3.c:1:
../../../libs/commons/string_utils.h:12:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from test3.h:13,
                 from test3.c:1:
BW_io.h:52:1: warning: "INSERTION" redefined
In file included from ../../../libs/containers/array_list.h:9,
                 from test3.h:9,
                 from test3.c:1:
../../../libs/commons/string_utils.h:11:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from test3.h:13,
                 from test3.c:1:
BW_io.h:162:1: warning: "tic" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from test3.h:12,
                 from test3.c:1:
../../../libs/commons/commons.h:67:1: warning: this is the location of the previous definition
In file included from bwt.h:13,
                 from test3.h:13,
                 from test3.c:1:
BW_io.h:167:1: warning: "toc" redefined
In file included from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from test3.h:12,
                 from test3.c:1:
../../../libs/commons/commons.h:72:1: warning: this is the location of the previous definition
test3.c: In function ‘main’:
test3.c:55: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
test3.c:55: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
test3.c:55: warning: format ‘%i’ expects type ‘int’, but argument 6 has type ‘size_t’
test3.c:55: warning: format ‘%i’ expects type ‘int’, but argument 7 has type ‘size_t’
