/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string_length’:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/alignment.c:487:5: warning: function returns address of local variable
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/alignment.c: In function ‘convert_to_quality_string’:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/alignment.c:498:5: warning: function returns address of local variable
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.c:1:0:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:8:14: warning: extra tokens at end of #ifdef directive
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.c:1:0:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:214:14: warning: extra tokens at end of #ifdef directive
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.c:199:14: warning: extra tokens at end of #ifdef directive
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.c: In function ‘print_bam1’:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.c:339:5: warning: format ‘%x’ expects type ‘unsigned int’, but argument 3 has type ‘struct bam1_t *’
array_list.c: In function ‘array_list_swap’:
array_list.c:379:16: warning: initialization makes integer from pointer without a cast
fastq_batch_reader.c: In function ‘fastq_batch_reader_single’:
fastq_batch_reader.c:61:3: warning: implicit declaration of function ‘omp_get_thread_num’
In file included from /usr/include/stdio.h:930:0,
                 from fastq_file.h:4,
                 from fastq_file.c:1:
/usr/include/bits/stdio2.h: In function ‘fastq_fread_index_positions’:
fastq_file.c:191:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:192:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:193:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
/usr/include/bits/stdio2.h: In function ‘fastq_fread_batch_max_size’:
fastq_file.c:120:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:121:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:122:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from /usr/include/string.h:642:0,
                 from fastq_file.h:5,
                 from fastq_file.c:1:
/usr/include/bits/string3.h: In function ‘fastq_fread_max_size’:
fastq_file.c:79:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:80:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:81:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c: In function ‘fastq_fread_num_reads’:
fastq_file.c:43:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:44:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:45:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from bwt.h:18:0,
                 from bwt.c:1:
BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
bwt.c: In function ‘read_cals_new’:
bwt.c:173:21: warning: assignment discards qualifiers from pointer target type
bwt.c: In function ‘read_cals_free’:
bwt.c:180:3: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
bwt.c: In function ‘bwt_generate_index_files’:
bwt.c:303:3: warning: implicit declaration of function ‘calculateBWT’
bwt.c:314:3: warning: implicit declaration of function ‘calculateC’
bwt.c:315:3: warning: implicit declaration of function ‘calculateO’
bwt.c:331:3: warning: implicit declaration of function ‘calculateR’
bwt.c:336:3: warning: implicit declaration of function ‘calculateSRcomp’
bwt.c: In function ‘bwt_map_inexact_seq’:
bwt.c:1220:8: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1226:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1229:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1232:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1239:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1242:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1245:3: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1252:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1255:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:1258:3: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:1263:8: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1263:8: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c: In function ‘print_item’:
bwt.c:1731:2: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:1731:2: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c: In function ‘bwt_generate_cal_list’:
bwt.c:2014:8: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:2014:8: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:2018:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:2018:3: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:2018:3: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:2018:3: warning: format ‘%d’ expects type ‘int’, but argument 6 has type ‘size_t’
bwt.c:2023:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:2023:3: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:2023:3: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:2023:3: warning: format ‘%d’ expects type ‘int’, but argument 6 has type ‘size_t’
bwt.c:2027:8: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:2027:8: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:2056:5: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:2056:5: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt.c:2056:5: warning: format ‘%d’ expects type ‘int’, but argument 4 has type ‘size_t’
bwt.c:2056:5: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
bwt.c:2066:5: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt.c:2066:5: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
BW_io.c: In function ‘load_exome_file’:
BW_io.c:685:10: warning: initialization makes integer from pointer without a cast
In file included from bwt.h:18:0,
                 from test-preprocess.c:1:
BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from bwt.h:18:0,
                 from test-search-cpu.c:1:
BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/../../../common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
test-search-cpu.c: In function ‘main’:
test-search-cpu.c:28:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
test-search-cpu.c:30:5: warning: format ‘%i’ expects type ‘int’, but argument 2 has type ‘size_t’
hpg_sw_bench.c: In function ‘run_emboss’:
hpg_sw_bench.c:273:3: warning: implicit declaration of function ‘omp_set_num_threads’
hpg_sw_bench.c:309:2: warning: implicit declaration of function ‘omp_get_thread_num’
fastq_batch_reader.c: In function ‘fastq_batch_reader_single’:
fastq_batch_reader.c:61:3: warning: implicit declaration of function ‘omp_get_thread_num’
In file included from /usr/include/stdio.h:930:0,
                 from fastq_file.h:4,
                 from fastq_file.c:1:
/usr/include/bits/stdio2.h: In function ‘fastq_fread_index_positions’:
fastq_file.c:191:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:192:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:193:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
/usr/include/bits/stdio2.h: In function ‘fastq_fread_batch_max_size’:
fastq_file.c:120:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:121:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:122:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from /usr/include/string.h:642:0,
                 from fastq_file.h:5,
                 from fastq_file.c:1:
/usr/include/bits/string3.h: In function ‘fastq_fread_max_size’:
fastq_file.c:79:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:80:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:81:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c: In function ‘fastq_fread_num_reads’:
fastq_file.c:43:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:44:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
fastq_file.c:45:14: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from main.c:17:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from batch_writer.h:17:0,
                 from main.c:19:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:8:14: warning: extra tokens at end of #ifdef directive
In file included from batch_writer.h:17:0,
                 from main.c:19:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:214:14: warning: extra tokens at end of #ifdef directive
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/sse.h:6:0,
                 from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/smith_waterman.h:8,
                 from sw_server.h:13,
                 from batch_writer.h:21,
                 from main.c:19:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:12:0: warning: "MAX" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:21:0: note: this is the location of the previous definition
main.c: In function ‘main’:
main.c:223:5: warning: implicit declaration of function ‘run_rna_aligner’
main.c:227:5: warning: implicit declaration of function ‘run_dna_aligner’
main.c: At top level:
main.c:270:6: warning: conflicting types for ‘run_dna_aligner’
main.c:227:5: note: previous implicit declaration of ‘run_dna_aligner’ was here
main.c: In function ‘run_dna_aligner’:
main.c:363:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
main.c:374:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
main.c:385:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
main.c:396:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
main.c: At top level:
main.c:414:6: warning: conflicting types for ‘run_rna_aligner’
main.c:223:5: note: previous implicit declaration of ‘run_rna_aligner’ was here
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/sse.h:6:0,
                 from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/smith_waterman.h:8,
                 from sw.h:15,
                 from sw.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:22:0: warning: "FLT_MAX" redefined
/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5.2/include/float.h:104:0: note: this is the location of the previous definition
sw.c: In function ‘smith_waterman_simd2’:
sw.c:296:25: warning: initialization from incompatible pointer type
sw.c:297:25: warning: initialization from incompatible pointer type
sw.c:431:22: warning: assignment from incompatible pointer type
sw.c:532:14: warning: assignment from incompatible pointer type
sw.c: In function ‘smith_waterman_simd’:
sw.c:651:3: warning: implicit declaration of function ‘reallocate_memory’
sw.c: In function ‘saveCharMatrixByIndex’:
sw.c:1208:5: warning: assignment from incompatible pointer type
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from sw_server.h:16,
                 from rna_server.h:6,
                 from rna_server.c:3:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from buffers.h:8:0,
                 from sw_server.h:16,
                 from rna_server.h:6,
                 from rna_server.c:3:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:21:0: warning: "MAX" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:12:0: note: this is the location of the previous definition
rna_server.c: In function ‘rna_server_omp_smith_waterman’:
rna_server.c:1262:3: warning: format ‘%i’ expects type ‘int’, but argument 2 has type ‘size_t’
rna_server.c:1262:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
rna_server.c:1262:3: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
rna_server.c:1262:3: warning: format ‘%i’ expects type ‘int’, but argument 5 has type ‘size_t’
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from rna_splice.h:8,
                 from rna_splice.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
rna_splice.c: In function ‘allocate_new_splice’:
rna_splice.c:108:3: warning: passing argument 2 of ‘cp_avltree_get’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/avl.h:116:7: note: expected ‘void *’ but argument is of type ‘unsigned int’
rna_splice.c:111:5: warning: passing argument 2 of ‘cp_avltree_insert’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/avl.h:112:7: note: expected ‘void *’ but argument is of type ‘unsigned int’
rna_splice.c:111:5: warning: passing argument 3 of ‘cp_avltree_insert’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/avl.h:112:7: note: expected ‘void *’ but argument is of type ‘unsigned int’
rna_splice.c: In function ‘process_avlnode_ends_in_order’:
rna_splice.c:189:209: warning: dereferencing ‘void *’ pointer
rna_splice.c:191:232: warning: dereferencing ‘void *’ pointer
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from bwt_server.h:7,
                 from bwt_server_cpu.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
bwt_server_cpu.c: In function ‘bwt_server_cpu’:
bwt_server_cpu.c:33:2: warning: implicit declaration of function ‘bwt_map_inexact_batch’
bwt_server_cpu.c:90:5: warning: format ‘%i’ expects type ‘int’, but argument 2 has type ‘size_t’
bwt_server_cpu.c:90:5: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
bwt_server_cpu.c:90:5: warning: format ‘%i’ expects type ‘int’, but argument 4 has type ‘size_t’
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from region_seeker.h:8,
                 from region_seeker.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
region_seeker.c: In function ‘region_seeker_server’:
region_seeker.c:56:7: warning: implicit declaration of function ‘bwt_map_exact_seeds_seq’
region_seeker.c: In function ‘apply_seeding’:
region_seeker.c:116:10: warning: assignment from incompatible pointer type
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from cal_seeker.h:4,
                 from cal_seeker.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
cal_seeker.c: In function ‘cal_seeker_server’:
cal_seeker.c:60:7: warning: implicit declaration of function ‘bwt_generate_cal_list_linkedlist’
cal_seeker.c:106:2: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
cal_seeker.c:155:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
cal_seeker.c: In function ‘apply_caling’:
cal_seeker.c:217:7: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct region_t *)’
cal_seeker.c:224:7: warning: passing argument 2 of ‘array_list_clear’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:77:5: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct region_t *)’
cal_seeker.c:225:7: warning: passing argument 2 of ‘array_list_clear’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:77:5: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from pair_server.h:13,
                 from pair_server.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
pair_server.c: In function ‘pair_server’:
pair_server.c:71:9: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
pair_server.c:76:4: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
pair_server.c:76:4: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
pair_server.c:76:4: warning: format ‘%d’ expects type ‘int’, but argument 5 has type ‘size_t’
pair_server.c:76:4: warning: format ‘%d’ expects type ‘int’, but argument 6 has type ‘size_t’
pair_server.c:136:3: warning: format ‘%d’ expects type ‘int’, but argument 2 has type ‘size_t’
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from sw_server.h:16,
                 from sw_server.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from buffers.h:8:0,
                 from sw_server.h:16,
                 from sw_server.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:21:0: warning: "MAX" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:12:0: note: this is the location of the previous definition
sw_server.c: In function ‘apply_sw’:
sw_server.c:474:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
In file included from batch_writer.h:17:0,
                 from batch_writer.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:8:14: warning: extra tokens at end of #ifdef directive
In file included from batch_writer.h:17:0,
                 from batch_writer.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/bam_file.h:214:14: warning: extra tokens at end of #ifdef directive
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from batch_writer.h:19,
                 from batch_writer.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/sse.h:6:0,
                 from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/smith_waterman.h:8,
                 from sw_server.h:13,
                 from batch_writer.h:21,
                 from batch_writer.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:12:0: warning: "MAX" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:21:0: note: this is the location of the previous definition
batch_writer.c: In function ‘batch_writer2’:
batch_writer.c:239:2: warning: format ‘%d’ expects type ‘int’, but argument 3 has type ‘size_t’
batch_writer.c:323:6: warning: passing argument 6 of ‘generate_cigar_str’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/bioformats/bam-sam/alignment.h:287:7: note: expected ‘short int *’ but argument is of type ‘size_t *’
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from batch_aligner.h:13,
                 from batch_aligner.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/sse.h:6:0,
                 from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/smith_waterman.h:8,
                 from sw_server.h:13,
                 from batch_aligner.h:19,
                 from batch_aligner.c:1:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/sw/macros.h:12:0: warning: "MAX" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:21:0: note: this is the location of the previous definition
In file included from /home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/bwt.h:18:0,
                 from buffers.h:8,
                 from buffers.c:7:
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:202:0: warning: "tic" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:67:0: note: this is the location of the previous definition
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:207:0: warning: "toc" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/commons/commons.h:72:0: note: this is the location of the previous definition
buffers.c:9:0: warning: "MAXLINE" redefined
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/bioinfo-libs/aligners/bwt/BW_io.h:43:0: note: this is the location of the previous definition
buffers.c: In function ‘region_batch_free’:
buffers.c:20:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct region_t *)’
buffers.c: In function ‘sw_batch_free’:
buffers.c:44:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct cal_t *)’
buffers.c: In function ‘cal_batch_free’:
buffers.c:121:5: warning: passing argument 2 of ‘array_list_free’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/common-libs/containers/array_list.h:65:6: note: expected ‘void (*)(void *)’ but argument is of type ‘void (*)(struct region_t *)’
timing.c: In function ‘timing_new’:
timing.c:27:3: warning: implicit declaration of function ‘calloc’
timing.c:27:31: warning: incompatible implicit declaration of built-in function ‘calloc’
timing.c:41:5: warning: implicit declaration of function ‘strdup’
timing.c:41:30: warning: assignment makes pointer from integer without a cast
timing.c: In function ‘timing_free’:
timing.c:59:5: warning: implicit declaration of function ‘free’
timing.c:59:45: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:61:43: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:62:41: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:63:41: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:66:40: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:68:38: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:69:37: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:70:36: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:72:37: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:73:39: warning: incompatible implicit declaration of built-in function ‘free’
timing.c:75:3: warning: incompatible implicit declaration of built-in function ‘free’
timing.c: In function ‘timing_display’:
timing.c:99:3: warning: implicit declaration of function ‘printf’
timing.c:99:3: warning: incompatible implicit declaration of built-in function ‘printf’
genome.c:12:35: warning: initializer element is not a constant expression
genome.c: In function ‘genome_new’:
genome.c:50:5: warning: passing argument 1 of ‘strlen’ from incompatible pointer type
/usr/include/string.h:399:15: note: expected ‘const char *’ but argument is of type ‘unsigned int *’
genome.c: In function ‘genome_get_chr_name’:
genome.c:106:3: warning: return from incompatible pointer type
genome.c: In function ‘load_hasthable_codes’:
genome.c:114:6: warning: passing argument 3 of ‘cp_hashtable_create’ from incompatible pointer type
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/hashtable.h:283:2: note: expected ‘cp_compare_fn’ but argument is of type ‘int (*)(const char *, const char *)’
genome.c:131:2: warning: passing argument 3 of ‘cp_hashtable_put’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/hashtable.h:477:7: note: expected ‘void *’ but argument is of type ‘unsigned char’
genome.c:144:7: warning: passing argument 3 of ‘cp_hashtable_put’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/hashtable.h:477:7: note: expected ‘void *’ but argument is of type ‘unsigned char’
genome.c:154:5: warning: passing argument 3 of ‘cp_hashtable_put’ makes pointer from integer without a cast
/home/hmartinez/appl/bioinfo-c/hpg-aligner/lib/cprops/hashtable.h:477:7: note: expected ‘void *’ but argument is of type ‘unsigned char’
genome.c: In function ‘code_binary_file_generator’:
genome.c:270:10: warning: cast from pointer to integer of different size
genome.c:285:13: warning: cast from pointer to integer of different size
genome.c: In function ‘load_binary_dna’:
genome.c:304:3: warning: format ‘%i’ expects type ‘int’, but argument 2 has type ‘size_t’
genome.c: In function ‘generate_codes’:
genome.c:476:3: warning: format ‘%i’ expects type ‘int’, but argument 2 has type ‘size_t’
genome.c:476:3: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
genome.c: In function ‘load_binary_dna’:
genome.c:307:8: warning: ignoring return value of ‘fread’, declared with attribute warn_unused_result
In file included from /usr/include/stdio.h:930:0,
                 from genome.h:4,
                 from genome.c:4:
/usr/include/bits/stdio2.h: In function ‘code_binary_file_generator’:
genome.c:264:10: warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result
statistics.c: In function ‘statistics_new’:
statistics.c:7:3: warning: implicit declaration of function ‘calloc’
statistics.c:7:38: warning: incompatible implicit declaration of built-in function ‘calloc’
statistics.c:20:5: warning: implicit declaration of function ‘strdup’
statistics.c:20:30: warning: assignment makes pointer from integer without a cast
statistics.c:26:34: warning: assignment makes pointer from integer without a cast
statistics.c: In function ‘statistics_free’:
statistics.c:35:7: warning: implicit declaration of function ‘free’
statistics.c:35:56: warning: incompatible implicit declaration of built-in function ‘free’
statistics.c:38:59: warning: incompatible implicit declaration of built-in function ‘free’
statistics.c: In function ‘statistics_display’:
statistics.c:73:2: warning: format ‘%-7d’ expects type ‘int’, but argument 3 has type ‘size_t’
statistics.c: In function ‘timing_and_statistics_display’:
statistics.c:106:7: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
statistics.c:109:7: warning: format ‘%i’ expects type ‘int’, but argument 3 has type ‘size_t’
