#include "hpg_sw.h"

double *sse_matrix_t, *sse_tracking_t;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void test();

// hpg-sw -q query_filename -r ref_filename -o output_filename -p gap_open_penalty -e gap_extend_penalty -s substitution_score_matrix -n number_of_threads -b number_of_reads_per_batch

int main(int argc, char *argv[]) {

    int c;
    char *q_filename = NULL, *r_filename = NULL;
    char *out_dirname = NULL, *out_filename = NULL, *subst_matrix_filename = NULL;
    int num_threads = 4, batch_size = 2048;
    float gap_open = 10.0f;
    float gap_extend = 0.5f;
    float match = 5.0f;
    float mismatch = -4.0f;

    while (1) {;
        static struct option long_options[] =
                {
                        {"query-file",               required_argument, 0, 'q'},
                        {"ref-file",                 required_argument, 0, 'r'},
                        {"output-dir",               optional_argument, 0, 'd'},
                        {"output-file",              optional_argument, 0, 'o'},
                        {"gap-open-penalty",         required_argument, 0, 'p'},
                        {"gap-extend-penalty",       required_argument, 0, 'e'},
                        {"match-score",              required_argument, 0, 'm'},
                        {"mismatch-penalty",         required_argument, 0, 'x'},
                        {"substitution-matrix-file", required_argument, 0, 's'},
                        {"num-threads",              required_argument, 0, 'n'},
                        {"reads-per-batch",          required_argument, 0, 'b'},
                        {"help",                     no_argument,       0, 'h'},
                        {0, 0, 0, 0}
                };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "q:r:d:o:p:e:m:x:s:n:b:h",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'q':
                q_filename = optarg;
                break;

            case 'r':
                r_filename = optarg;
                break;

            case 'd':
                out_dirname = optarg;
                break;

            case 'o':
                out_filename = optarg;
                break;

            case 'p':
                sscanf(optarg, "%f", &gap_open);
                break;

            case 'e':
                sscanf(optarg, "%f", &gap_extend);
                break;

            case 'm':
                sscanf(optarg, "%f", &match);
                break;

            case 'x':
                sscanf(optarg, "%f", &mismatch);
                break;

            case 's':
                subst_matrix_filename = optarg;
                break;

            case 'n':
                sscanf(optarg, "%i", &num_threads);
                break;

            case 'b':
                sscanf(optarg, "%i", &batch_size);
                break;

            case '?':
            case 'h':
            default:
                display_usage(NULL);
        }
    }

    // sanity cheking
    if (q_filename == NULL) {
        display_usage("\nERROR: Query file name is missing (option --query-file or -q)\n");
    } else if (!file_exists(q_filename)) {
        printf("\nERROR: Query file %s does not exist\n", q_filename);
        exit(-1);
    }

    printf("query-file = %s\n", q_filename);

    if (r_filename == NULL) {
        display_usage("\nERROR: Reference file name is missing (option --ref-file or -r)\n");
    } else if (!file_exists(r_filename)) {
        printf("\nERROR: Reference file %s does not exist\n", r_filename);
        exit(-1);
    }
    printf("ref-file = %s\n", r_filename);

    if (out_dirname == NULL) {
      //display_usage("\nERROR: Output directory name is missing (option --output-dir or -d)\n");
    } else if (!file_exists(out_dirname)) {
        printf("\nERROR: Output directory %s does not exist\n", out_dirname);
        exit(-1);
    }
    //printf("output-dir = %s\n", out_dirname);

    //    if (out_filename == NULL) {
    //        display_usage("\nERROR: Ouptut file name is missing (option --ouput-file or -o)\n");
    //    }
    //printf("output-file = %s\n", out_filename);

    printf("match-penalty = %f\n", match);
    printf("mismatch-penalty = %f\n", mismatch);

    if (gap_open < 0.0 || gap_open > 100.0) {
        display_usage("\nERROR: Penalty for the gap openning is not valid, it must be from 0.0 to 100.0 (option --gap-open-penalty or -p)\n");
    }
    printf("gap-open-penalty = %f\n", gap_open);

    if (gap_extend < 0.0 || gap_extend > 10.0) {
        display_usage("\nERROR: Penalty for the gap extending is not valid, it must be from 0.0 to 10.0 (option --gap-extend-penalty or -e)\n");
    }
    printf("gap-extend-penalty = %f\n", gap_extend);


    if (subst_matrix_filename == NULL) {
        display_usage("\nERROR: Substitution score matrix file name is missing (option --substitution-matrix-filename or -s)\n");
    } else if (!file_exists(subst_matrix_filename)) {
        printf("\nERROR: Substitution score matrix file %s does not exist\n", subst_matrix_filename);
        exit(-1);
    }
    printf("substitution-matrix-filename = %s\n", subst_matrix_filename);

    if (num_threads <= 0) {
        display_usage("\nERROR: Number of threads is not valid, it must be greater than 0 (option --num-threads or -n)\n");
    }
    printf("num-threads = %i\n", num_threads);

    if (batch_size <= 0) {
        display_usage("\nERROR: Number of reads per batch is not valid, it must be greater than 0 (option --reads-per-batch or -b)\n");
    }
    printf("reads-per-batch = %i\n", batch_size);

    char *sse_out_filename = NULL;
    if (out_filename == NULL) {
      printf("ATENTION: Missing output filename,runnig Smith-Waterman and backtracking for timing purposes\n");
    } else {
      sse_out_filename = (char *) malloc(2048 * sizeof(char));
      sprintf(sse_out_filename, "%s/%s", out_dirname, out_filename);
      printf("out filename (full path) = %s\n", sse_out_filename);
    }

    sse_matrix_t = (double *) calloc(num_threads, sizeof(double));
    sse_tracking_t = (double *) calloc(num_threads, sizeof(double));

    // run_sw (was the original call)
    run_sw_omp(q_filename, r_filename,
	       match, mismatch, gap_open, gap_extend, subst_matrix_filename, batch_size,
	       num_threads, sse_out_filename);

    if (out_filename) {
      free(sse_out_filename);
    }
    free(sse_matrix_t);
    free(sse_tracking_t);

    printf("Done.\n");
}

//-------------------------------------------------------------------------
void run_sw(char *q_filename, char *r_filename,
	    float match, float mismatch, float gap_open, float gap_extend,
	    char *matrix_filename, int batch_size, int num_threads,
	    char *out_filename) {

  printf("------------------------------------------\n");
  printf("Running function %s\n", __FUNCTION__);
  printf("------------------------------------------\n");

    const int max_length = 2048;

    double start = 0.0f, elapsed = 0.0f;
    double start_read = 0.0f, read = 0.0f;
    double start_write = 0.0f, write = 0.0f;
    double start_memory = 0.0f, memory = 0.0f;
    double sse_t = 0.0f, partial_t = 0.0f;

    sw_optarg_t *optarg_p = sw_optarg_new(match, mismatch, gap_open, gap_extend, matrix_filename);
    sw_multi_output_t *output_p;

    FILE *q_file = fopen(q_filename, "r");
    FILE *r_file = fopen(r_filename, "r");
    FILE *out_file = fopen(out_filename, "w");

    start = sw_tic();
#ifdef TIMING
    start_memory = start;
#endif
    char *q[batch_size], *r[batch_size];
    for (int i = 0; i < batch_size; i++) {
        q[i] = (char *) calloc(max_length, sizeof(char));
        r[i] = (char *) calloc(max_length, sizeof(char));
    }
#ifdef TIMING
    memory += sw_toc(start_memory);
#endif

    int count = 0, batches = 0, num_queries;

    while (1) {
        num_queries = 0;

#ifdef TIMING
	start_read = sw_tic();
#endif
        // read queries
        for (int i = 0; i < batch_size; i++) {
            if (fgets(q[i], max_length, q_file) == NULL) { break; }
            str_trim(q[i]);
            num_queries++;
            count++;
        }

        // exit if no queries
        if (num_queries == 0) break;


        // read references
        for (int i = 0; i < num_queries; i++) {
            if (fgets(r[i], max_length, r_file) == NULL) { break; }
            str_trim(r[i]);
        }
#ifdef TIMING
	read += sw_toc(start_read);
#endif

#ifdef TIMING
	start_memory = sw_tic();
#endif
        output_p = sw_multi_output_new(num_queries);
#ifdef TIMING
	memory += sw_toc(start_memory);
#endif

        // call smith-waterman
        partial_t = sw_tic();
        smith_waterman_mqmr(q, r, num_queries, optarg_p, num_threads, output_p);
        sse_t += sw_toc(partial_t);

        // save results
#ifdef TIMING
	start_write = sw_tic();
#endif
        sw_multi_output_save(num_queries, output_p, out_file);
#ifdef TIMING
	write += sw_toc(start_write);
#endif

#ifdef TIMING
	start_memory = sw_tic();
#endif
        sw_multi_output_free(output_p);
#ifdef TIMING
	memory += sw_toc(start_memory);
#endif

        batches++;
        //printf("%i batches\n", batches);
        //    if (batches == 7) break;
        //break;
    }

#ifdef TIMING
    start_memory = sw_tic();
#endif
    // free memory
    sw_optarg_free(optarg_p);
    for (int i = 0; i < batch_size; i++) {
        free(q[i]);
        free(r[i]);
    }
#ifdef TIMING
    memory += sw_toc(start_memory);
#endif

    // close files
    fclose(q_file);
    fclose(r_file);
    fclose(out_file);

    elapsed += sw_toc(start);

#ifdef TIMING
    double max_sse = 0.0f;
    for (int i = 0; i < num_threads; i++) {
        if (sse_matrix_t[i] + sse_tracking_t[i] > max_sse)
            max_sse = sse_matrix_t[i] + sse_tracking_t[i];
    }
#endif

#if defined (SW_AVX512) || defined (__AVX512__)
    printf("AVX512 (time in seconds)\n");
#else
#if defined (SW_AVX2) || defined (__AVX2__)
  printf("AVX2 (time in seconds)\n");
#else
  printf("SSE (time in seconds)\n");
#endif
#endif

#ifdef TIMING
    printf("Thread: <matrix creation time> <backtracing time>\n");
    for(int i = 0; i < num_threads ; i++) {
      printf("\tThread %i:\t%0.3fs\t%0.3fs\n", i, sse_matrix_t[i], sse_tracking_t[i]);
    }
    printf("\tMax. time: %0.3fs\n", max_sse);
#endif
#ifdef TIMING
    printf("Memory mng. time    : %0.3f s\n", memory);
    printf("Read sequences time : %0.3f s\n", read);
    printf("Write results time  : %0.3f s\n", write);
#endif
    printf("Alignment time      : %0.3f s (2 x %i seqs; %i threads)\n", sse_t, count, num_threads);
    printf("\n");
    printf("Total hpg-sw time   : %0.3fs\n", elapsed);
}


//-------------------------------------------------------------------------

void display_usage(char *msg) {
    if (msg != NULL) printf("%s", msg);
    printf("Usage: ./hpg-sw -q query_filename -r ref_filename -d output_dirname -o output_filename -p gap_open_penalty -e gap_extend_penalty -s substitution_score_matrix -n number_of_threads -b number_of_reads_per_batch\n");

    printf("\t--query-file, -q, input file name containing the queries to align\n");
    printf("\t--ref-file, -r, input file name containing the references to align\n");
    printf("\t--output-dir, -d, output directory where the results will be saved\n");
    printf("\t--output-file, -o, output file name where the alignments will be saved\n");
    printf("\t--match-score, -m, score for matches; default 5.0\n");
    printf("\t--mismatch-penalty, -x, penalty for mismatches; default -4.0\n");
    printf("\t--gap-open-penalty, -p, penalty for the gap openning: from 0.0 to 100.0\n");
    printf("\t--gap-extend-penalty, -e, penalty for the gap extending: from 0.0 to 10.0\n");
    printf("\t--substitution-matrix, -s, substitution score matrix name, for DNA: dnafull, for proteins: blosum50, blosum62, blosum80\n");
    printf("\t--num-threads, -n, number of threads\n");
    printf("\t--reads-per-batch, -b, number of reads per batch\n");

    exit(-1);
}

//-------------------------------------------------------------------------

int file_exists(const char *filename) {
    FILE *file = NULL;
    if (file = fopen(filename, "r")) {
        fclose(file);
        return 1;
    }
    return 0;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
