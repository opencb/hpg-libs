In file included from ../../../libs/containers/array_list.h:9:0,
                 from bwt.h:9,
                 from bwt.c:1:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../hpg-aligner/samtools-0.1.18/bam.h:52:0,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from bwt.h:12,
                 from bwt.c:1:
../../../hpg-aligner/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
../../../hpg-aligner/samtools-0.1.18/bgzf.h:149:9: warning: implicit declaration of function ‘ftello’
In file included from bwt.h:13:0,
                 from bwt.c:1:
BW_io.h: At top level:
BW_io.h:50:0: warning: "DELETION" redefined
../../../libs/commons/string_utils.h:12:0: note: this is the location of the previous definition
BW_io.h:52:0: warning: "INSERTION" redefined
../../../libs/commons/string_utils.h:11:0: note: this is the location of the previous definition
BW_io.h:162:0: warning: "tic" redefined
../../../libs/commons/commons.h:67:0: note: this is the location of the previous definition
BW_io.h:167:0: warning: "toc" redefined
../../../libs/commons/commons.h:72:0: note: this is the location of the previous definition
bwt.c: In function ‘read_cals_new’:
bwt.c:88:18: warning: assignment discards qualifiers from pointer target type
bwt.c: In function ‘read_cals_free’:
bwt.c:93:3: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
bwt.c: In function ‘bwt_index_new’:
bwt.c:144:3: warning: implicit declaration of function ‘strdup’
bwt.c:144:18: warning: assignment makes pointer from integer without a cast
bwt.c: In function ‘bwt_map_exact_batch_cpu’:
bwt.c:527:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct fastq_read_t *)’
bwt.c: In function ‘bwt_map_inexact_seqs_by_pos_cpu’:
bwt.c:968:7: warning: passing argument 1 of ‘bwt_map_inexact_seq_cpu’ makes pointer from integer without a cast
bwt.c:742:14: note: expected ‘char *’ but argument is of type ‘char’
bwt.c: In function ‘bwt_map_inexact_batch_cpu’:
bwt.c:1159:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct fastq_read_t *)’
bwt.c: In function ‘bwt_map_inexact_seeds_seq_cpu’:
bwt.c:1185:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c: In function ‘create_seeds_from_seqs_cpu’:
bwt.c:1310:11: warning: passing argument 4 of ‘create_seeds_from_seq_cpu’ makes pointer from integer without a cast
bwt.c:1219:9: note: expected ‘unsigned int *’ but argument is of type ‘unsigned int’
bwt.c: In function ‘create_seeds_from_batch_cpu’:
bwt.c:1341:13: warning: passing argument 4 of ‘create_seeds_from_seq_cpu’ makes pointer from integer without a cast
bwt.c:1219:9: note: expected ‘unsigned int *’ but argument is of type ‘unsigned int’
bwt.c: In function ‘bwt_generate_cal_list_cpu’:
bwt.c:1382:5: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1382:5: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1404:8: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1404:8: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1404:8: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1408:8: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1408:8: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1408:8: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1411:6: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1411:6: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1439:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1439:3: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1439:3: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1439:3: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:1448:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1448:3: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1480:5: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1480:5: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1480:5: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1484:7: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1484:7: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c: In function ‘bwt_find_cals_from_seqs_cpu’:
bwt.c:1610:7: warning: passing argument 1 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65:6: note: expected ‘struct array_list_t *’ but argument is of type ‘struct array_list_t **’
bwt.c:1617:9: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
../../../libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct region_t *)’
BW_io.c: In function ‘load_exome_file’:
BW_io.c:651:10: warning: initialization makes integer from pointer without a cast
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/commons/file_utils.h:13,
                 from ../../../libs/bioformats/fastq/fastq_file.h:8,
                 from ../../../libs/bioformats/fastq/fastq_file.c:7:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from /usr/include/stdio.h:930:0,
                 from ../../../libs/bioformats/fastq/fastq_file.c:2:
/usr/include/bits/stdio2.h: In function ‘fastq_fread_index_positions’:
../../../libs/bioformats/fastq/fastq_file.c:201:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:202:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:203:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
/usr/include/bits/stdio2.h: In function ‘fastq_fread_batch_max_size’:
../../../libs/bioformats/fastq/fastq_file.c:130:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:131:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:132:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from /usr/include/string.h:642:0,
                 from ../../../libs/bioformats/fastq/fastq_file.c:3:
/usr/include/bits/string3.h: In function ‘fastq_fread_max_size’:
../../../libs/bioformats/fastq/fastq_file.c:89:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:90:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:91:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c: In function ‘fastq_fread_num_reads’:
../../../libs/bioformats/fastq/fastq_file.c:53:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:54:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/bioformats/fastq/fastq_file.c:55:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from ../../../libs/bioformats/fastq/fastq_batch.c:7:0:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.h:11,
                 from ../../../libs/bioformats/fastq/fastq_batch_reader.c:11:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../hpg-aligner/samtools-0.1.18/bam.h:52:0,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
../../../hpg-aligner/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
../../../hpg-aligner/samtools-0.1.18/bgzf.h:149:9: warning: implicit declaration of function ‘ftello’
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/bioformats/bam-sam/alignment.h:9,
                 from ../../../libs/bioformats/bam-sam/alignment.c:2:
../../../libs/commons/string_utils.h: At top level:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../libs/bioformats/bam-sam/alignment.c:3:0:
../../../libs/bioformats/bam-sam/bam_file.h:7:14: warning: extra tokens at end of #ifdef directive
In file included from ../../../libs/bioformats/bam-sam/alignment.c:3:0:
../../../libs/bioformats/bam-sam/bam_file.h:211:14: warning: extra tokens at end of #ifdef directive
../../../libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string_length’:
../../../libs/bioformats/bam-sam/alignment.c:416:5: warning: function returns address of local variable
../../../libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string’:
../../../libs/bioformats/bam-sam/alignment.c:427:5: warning: function returns address of local variable
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/commons/log.c:1:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/commons/commons.h:8,
                 from ../../../libs/commons/system_utils.h:10,
                 from ../../../libs/commons/system_utils.c:4:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
../../../libs/commons/system_utils.c: In function ‘get_optimal_cpu_num_threads’:
../../../libs/commons/system_utils.c:69:5: warning: implicit declaration of function ‘popen’
../../../libs/commons/system_utils.c:69:22: warning: assignment makes pointer from integer without a cast
../../../libs/commons/system_utils.c:82:5: warning: implicit declaration of function ‘pclose’
../../../libs/commons/system_utils.c:75:10: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
../../../libs/commons/system_utils.c: In function ‘get_free_memory’:
../../../libs/commons/system_utils.c:16:11: warning: ignoring return value of ‘system’, declared with attribute warn_unused_result
../../../libs/commons/system_utils.c:23:10: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from ../../../libs/commons/string_utils.c:1:0:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
../../../libs/commons/string_utils.c: In function ‘splitn’:
../../../libs/commons/string_utils.c:342:5: warning: implicit declaration of function ‘strsep’
../../../libs/commons/string_utils.c:342:19: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c:343:9: warning: implicit declaration of function ‘strdup’
../../../libs/commons/string_utils.c:343:23: warning: assignment makes pointer from integer without a cast
../../../libs/commons/string_utils.c: In function ‘get_to_first_blank’:
../../../libs/commons/string_utils.c:441:5: warning: return makes integer from pointer without a cast
In file included from ../../../libs/containers/array_list.h:9:0,
                 from ../../../libs/containers/array_list.c:4:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../libs/commons/log.h:10:0,
                 from ../../../libs/containers/list.h:8,
                 from ../../../libs/containers/list.c:5:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../libs/containers/array_list.h:9:0,
                 from test3.h:9,
                 from test3.c:1:
../../../libs/commons/string_utils.h:17:0: warning: "PADDING" redefined
../../../libs/commons/string_utils.h:15:0: note: this is the location of the previous definition
In file included from ../../../hpg-aligner/samtools-0.1.18/bam.h:52:0,
                 from ../../../libs/bioformats/bam-sam/alignment.h:7,
                 from test3.h:12,
                 from test3.c:1:
../../../hpg-aligner/samtools-0.1.18/bgzf.h: In function ‘bgzf_getc’:
../../../hpg-aligner/samtools-0.1.18/bgzf.h:149:9: warning: implicit declaration of function ‘ftello’
In file included from bwt.h:13:0,
                 from test3.h:13,
                 from test3.c:1:
BW_io.h: At top level:
BW_io.h:50:0: warning: "DELETION" redefined
../../../libs/commons/string_utils.h:12:0: note: this is the location of the previous definition
BW_io.h:52:0: warning: "INSERTION" redefined
../../../libs/commons/string_utils.h:11:0: note: this is the location of the previous definition
BW_io.h:162:0: warning: "tic" redefined
../../../libs/commons/commons.h:67:0: note: this is the location of the previous definition
BW_io.h:167:0: warning: "toc" redefined
../../../libs/commons/commons.h:72:0: note: this is the location of the previous definition
test3.c: In function ‘main’:
test3.c:56:5: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
test3.c:56:5: warning: format ‘%i’ expects type ‘int’, but argument 6 has type ‘size_t’
test3.c:56:5: warning: format ‘%i’ expects type ‘int’, but argument 7 has type ‘size_t’
