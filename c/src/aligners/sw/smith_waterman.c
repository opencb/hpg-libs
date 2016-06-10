#include "smith_waterman.h"

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#ifdef TIMING
extern double *sse_matrix_t, *sse_tracking_t;
#endif // TIMING

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_mqmr(char **query_p, char **ref_p, unsigned int num_queries,
                         sw_optarg_t *optarg_p, unsigned int num_threads,
                         sw_multi_output_t *output_p) {
#ifdef __AVX2__
    const unsigned int simd_depth = 8;
#else
    const unsigned int simd_depth = 4;
#endif


    if (output_p == NULL) {
        printf("Error: output buffer is null.\n");
        exit(-1);
    }
    if (output_p->num_queries < num_queries) {
        printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n",
               num_queries, output_p->num_queries);
        exit(-1);
    }

    if (num_threads == 1) {
        int tid = omp_get_thread_num();

#ifdef TIMING
        double partial_t;
#endif // TIMING
        char *q_aux = NULL, *r_aux = NULL;
        int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
        char *q[simd_depth], *r[simd_depth];
        int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
        float *H = NULL, *F = NULL;
        int *C = NULL;

        float match = optarg_p->match;
        float mismatch = optarg_p->mismatch;
        float gap_open = optarg_p->gap_open;
        float gap_extend = optarg_p->gap_extend;
        float *score_p = output_p->score_p;

        //printf("num queries = %i\n", num_queries);
        //for(int k = 0; k < 8; k++) printf("%0.2f ", score_p[k]);
        //printf("\n");

        depth = 0;
        for (unsigned int i = 0; i < num_queries; i++) {
            //printf("smith_waterman.c: query: #%i\n", i);
            len = strlen(query_p[i]);
            if (len > max_q_len) max_q_len = len;
            q_lens[depth] = len;
            q[depth] = query_p[i];

            len = strlen(ref_p[i]);
            if (len > max_r_len) max_r_len = len;
            r_lens[depth] = len;
            r[depth] = ref_p[i];

            depth++;
            if (depth == simd_depth) {

                index = i - depth + 1;

                reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
                //printf("-----> max_q_len = %i, max_r_eln = %i, simd_depth = %i\n", max_q_len, max_r_len, simd_depth);
                //printf("-----> H_size = %i, F_size = %i, aux_size = %i\n", H_size, F_size, aux_size);
                //printf("num_queries = %i, gap_open = %0.2f, gap_extend = %0.2f\n", depth, gap_open, gap_extend);

                // generating score matrix
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
#ifdef SW_VERSION_1
                avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                avx2_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                             match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif // SW_VERSION_1
#else
#ifdef SW_VERSION_1
                sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                sse_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif // SW_VERSION_1
#endif // __AVX2__

                //printf("start avx2_matrix (index %i)\n", index);
                //printf("end avx2_matrix\n");

#ifdef TIMING
                sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING
                // tracebacking
                simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               gap_open, gap_extend, H, C, &score_p[index],
                               &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                               &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                               q_aux, r_aux);
#ifdef TIMING
                sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
                depth = 0;
                max_q_len = 0;
                max_r_len = 0;
            }
        }

        //    printf("depth = %i\n", depth);

        if (depth > 0) {

            //float max_score[simd_depth];
            float *max_score = (float *) _mm_malloc(simd_depth * sizeof(float), 32);

            for (unsigned int i = depth; i < simd_depth; i++) {
                q[i] = q[0];
                q_lens[i] = q_lens[0];

                r[i] = r[0];
                r_lens[i] = r_lens[0];
            }

            reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

            index = num_queries - depth;

            // generating score matrix
#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
#ifdef SW_VERSION_1
            avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                        optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
            avx2_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                         match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#else
#ifdef SW_VERSION_1
            sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                       optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
            sse_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                        match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#endif // __AVX2__

#ifdef TIMING
            sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING
            // tracebacking
            simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           gap_open, gap_extend, H, C, max_score,
                           &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                           &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                           q_aux, r_aux);
#ifdef TIMING
            sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

            for (unsigned int i = 0; i < depth; i++) {
                score_p[index +  i] = max_score[i];
            }
            _mm_free(max_score);

        }
        //    printf("Free 1\n");
        // free memory
        if (H != NULL) _mm_free(H);
        if (C != NULL) _mm_free(C);
        if (F != NULL) _mm_free(F);
        if (q_aux != NULL) free(q_aux);
        if (r_aux != NULL) free(r_aux);

    } else { // multi-thread

        unsigned int num_packs = num_queries / simd_depth;
        if (num_queries % simd_depth) {
            num_packs++;
        }
        unsigned int packs_per_thread = num_packs / num_threads;
        if (num_packs % num_threads) {
            packs_per_thread++;
        }

        //printf("num_packs = %i, packs_per_thread = %i\n", num_packs, packs_per_thread);

#pragma omp parallel num_threads(num_threads) shared(num_packs, packs_per_thread, optarg_p, output_p)
        {
#ifdef TIMING
            double partial_t = 0.0;
#endif // TIMING
            unsigned int tid = omp_get_thread_num();

            unsigned int first_index = tid * packs_per_thread * simd_depth;
            unsigned int last_index = first_index + (packs_per_thread * simd_depth);
            if (last_index > num_queries) last_index = num_queries;

            //printf("tid = %i, from %i to %i\n", tid, first_index, last_index);

            char *q_aux = NULL, *r_aux = NULL;
            int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
            char *q[simd_depth], *r[simd_depth];
            int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
            float *H = NULL, *F = NULL;
            int *C = NULL;

            float match = optarg_p->match;
            float mismatch = optarg_p->mismatch;
            float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
            float *score_p = output_p->score_p;

            depth = 0;
            for (unsigned int i = first_index; i < last_index; i++) {

                len = strlen(query_p[i]);
                if (len > max_q_len) max_q_len = len;
                q_lens[depth] = len;
                q[depth] = query_p[i];

                len = strlen(ref_p[i]);
                if (len > max_r_len) max_r_len = len;
                r_lens[depth] = len;
                r[depth] = ref_p[i];

                depth++;
                if (depth == simd_depth) {

                    index = i - depth + 1;

                    reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

                    // generating score matrix
#ifdef TIMING
                    partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
#ifdef SW_VERSION_1
                    avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                    avx2_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                 match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif // SW_VERSION_1
#else
#ifdef SW_VERSION_1
                    sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                    sse_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif // SW_VERSION_1
#endif // __AVX2__

#ifdef TIMING
                    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

                    // tracebacking
#ifdef TIMING
                    partial_t = sw_tic();
#endif // TIMING
                    simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                   gap_open, gap_extend, H, C, &score_p[index],
                                   &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                                   &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                                   q_aux, r_aux);
#ifdef TIMING
                    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

                    depth = 0;
                    max_q_len = 0;
                    max_r_len = 0;
                }
            }

            //	printf("depth = %i\n", depth);

            if (depth > 0) {

                //float max_score[simd_depth];
                float *max_score = (float *) _mm_malloc(simd_depth * sizeof(float), 32);

                for (unsigned int i = depth; i < simd_depth; i++) {
                    q[i] = q[0];
                    q_lens[i] = q_lens[0];

                    r[i] = r[0];
                    r_lens[i] = r_lens[0];
                }

                reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

                index = num_queries - depth;

                // generating score matrix
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
#ifdef SW_VERSION_1
                avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
                avx2_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                             match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#else
#ifdef SW_VERSION_1
                sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
                sse_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#endif // __AVX2__

#ifdef TIMING
                sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

                // tracebacking
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING
                simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               gap_open, gap_extend, H, C, max_score,
                               &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                               &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                               q_aux, r_aux);
#ifdef TIMING
                sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

                for (unsigned int i = 0; i < depth; i++) {
                    score_p[index +  i] = max_score[i];
                }
                _mm_free(max_score);
            }

            // free memory
            if (H != NULL) _mm_free(H);
            if (C != NULL) _mm_free(C);
            if (F != NULL) _mm_free(F);
            if (q_aux != NULL) free(q_aux);
            if (r_aux != NULL) free(r_aux);

        } // end #pragma omp parallel
    }
}

//------------------------------------------------------------------------------------

void smith_waterman_mqsr(char **query_p, char *ref_p, unsigned int num_queries,
                         sw_optarg_t *optarg_p, unsigned int num_threads,
                         sw_multi_output_t *output_p) {

#ifdef __AVX2__
    const unsigned int simd_depth = 8;
#else
    const unsigned int simd_depth = 4;
#endif

#ifdef TIMING
    double partial_t;
#endif // TIMING

    if (output_p == NULL) {
        printf("Error: output buffer is null.\n");
        exit(-1);
    }
    if (output_p->num_queries < num_queries) {
        printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n",
               num_queries, output_p->num_queries);
        exit(-1);
    }

    int ref_len = strlen(ref_p);

    if (num_threads == 1) {

        int tid = omp_get_thread_num();

        char *q_aux = NULL, *r_aux = NULL;
        int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = ref_len;
        char *q[simd_depth], *r[simd_depth];
        int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
        float *H = NULL, *F = NULL;
        int *C = NULL;

        float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
        float *score_p = output_p->score_p;

        depth = 0;
        for (unsigned int i = 0; i < num_queries; i++) {

            len = strlen(query_p[i]);
            if (len > max_q_len) max_q_len = len;
            q_lens[depth] = len;
            q[depth] = query_p[i];

            r_lens[depth] = ref_len;
            r[depth] = ref_p;

            depth++;
            if (depth == simd_depth) {

                index = i - depth + 1;

                reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
                //printf("-----> max_q_len = %i, max_r_eln = %i, simd_depth = %i\n", max_q_len, max_r_len, simd_depth);
                //printf("-----> H_size = %i, F_size = %i, aux_size = %i\n", H_size, F_size, aux_size);
                //printf("num_queries = %i, gap_open = %0.2f, gap_extend = %0.2f\n", depth, gap_open, gap_extend);

                // generating score matrix
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING
#ifdef __AVX2__
                avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif
#ifdef TIMING
                sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

                // tracebacking
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING
                simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               gap_open, gap_extend, H, C, &score_p[index],
                               &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                               &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                               q_aux, r_aux);
#ifdef TIMING
                sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

                depth = 0;
                max_q_len = 0;
                max_r_len = 0;
            }
        }

        //    printf("depth = %i\n", depth);

        if (depth > 0) {

            //float max_score[simd_depth];
            float *max_score = (float *) _mm_malloc(simd_depth * sizeof(float), 32);

            for (unsigned int i = depth; i < simd_depth; i++) {
                q[i] = q[0];
                q_lens[i] = q_lens[0];

                r[i] = r[0];
                r_lens[i] = r_lens[0];
            }

            reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

            index = num_queries - depth;

            // generating score matrix
#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
            avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                        optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
            sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                       optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#endif

#ifdef TIMING
            sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

            // tracebacking
#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING
            simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           gap_open, gap_extend, H, C, max_score,
                           &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                           &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                           q_aux, r_aux);
#ifdef TIMING
            sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

            for (unsigned int i = 0; i < depth; i++) {
                score_p[index +  i] = max_score[i];
            }
            _mm_free(max_score);
        }

        // free memory
        if (H != NULL) _mm_free(H);
        if (C != NULL) _mm_free(C);
        if (F != NULL) _mm_free(F);
        if (q_aux != NULL) free(q_aux);
        if (r_aux != NULL) free(r_aux);

    } else { // multi-thread

        unsigned int num_packs = num_queries / simd_depth;
        if (num_queries % simd_depth) {
            num_packs++;
        }
        unsigned int packs_per_thread = num_packs / num_threads;
        if (num_packs % num_threads) {
            packs_per_thread++;
        }

        //printf("num_packs = %i, packs_per_thread = %i\n", num_packs, packs_per_thread);

#pragma omp parallel num_threads(num_threads) shared(num_packs, packs_per_thread, optarg_p, output_p)
        {
            unsigned int tid = omp_get_thread_num();

            unsigned int first_index = tid * packs_per_thread * simd_depth;
            unsigned int last_index = first_index + (packs_per_thread * simd_depth);
            if (last_index > num_queries) last_index = num_queries;

            //printf("tid = %i, from %i to %i\n", tid, first_index, last_index);

            char *q_aux = NULL, *r_aux = NULL;
            int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = ref_len;
            char *q[simd_depth], *r[simd_depth];
            int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
            float *H = NULL, *F = NULL;
            int *C = NULL;

            float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
            float *score_p = output_p->score_p;

            depth = 0;
            for (unsigned int i = first_index; i < last_index; i++) {

                len = strlen(query_p[i]);
                if (len > max_q_len) max_q_len = len;
                q_lens[depth] = len;
                q[depth] = query_p[i];

                r_lens[depth] = ref_len;
                r[depth] = ref_p;

                depth++;
                if (depth == simd_depth) {

                    index = i - depth + 1;

                    reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

                    // generating score matrix
#ifdef TIMING
                    partial_t = sw_tic();
#endif // TIMING
#ifdef __AVX2__
                    avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
                    sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#endif
#ifdef TIMING
                    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

                    // tracebacking
#ifdef TIMING
                    partial_t = sw_tic();
#endif // TIMING
                    simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                                   gap_open, gap_extend, H, C, &score_p[index],
                                   &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                                   &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                                   q_aux, r_aux);
#ifdef TIMING
                    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

                    depth = 0;
                    max_q_len = 0;
                    max_r_len = 0;
                }
            }

            //	printf("depth = %i\n", depth);

            if (depth > 0) {

                //float max_score[simd_depth];
                float *max_score = (float *) _mm_malloc(simd_depth * sizeof(float), 32);

                for (unsigned int i = depth; i < simd_depth; i++) {
                    q[i] = q[0];
                    q_lens[i] = q_lens[0];

                    r[i] = r[0];
                    r_lens[i] = r_lens[0];
                }

                reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

                index = num_queries - depth;

                // generating score matrix
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
                avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                            optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
                sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#endif

#ifdef TIMING
                sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

                // tracebacking
#ifdef TIMING
                partial_t = sw_tic();
#endif // TIMING
                simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                               gap_open, gap_extend, H, C, max_score,
                               &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                               &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                               q_aux, r_aux);
#ifdef TIMING
                sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

                for (unsigned int i = 0; i < depth; i++) {
                    score_p[index +  i] = max_score[i];
                }
                _mm_free(max_score);
            }

            // free memory
            if (H != NULL) _mm_free(H);
            if (C != NULL) _mm_free(C);
            if (F != NULL) _mm_free(F);
            if (q_aux != NULL) free(q_aux);
            if (r_aux != NULL) free(r_aux);

        } // end #pragma omp parallel
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void reallocate_memory(int max_q_len, int max_r_len, int simd_depth,
                       int *H_size, float **H, int **C, int *F_size, float **F,
                       int *aux_size, char **q_aux, char **r_aux) {

    int size_h, size_c, size_f;
    unsigned int matrix_size = max_q_len * max_r_len;

    size_h = simd_depth * matrix_size * sizeof(float);
    size_c = simd_depth * matrix_size * sizeof(int);
    if (matrix_size > *H_size) {
        if (*H_size > 0) {
            _mm_free(*H);
            _mm_free(*C);
        }

        *H = (float *) _mm_malloc(size_h, 32);
        *C = (int *) _mm_malloc(size_c, 32);
        //printf("new H %x, C %x\n", *H, *C);
        *H_size = matrix_size;
    }
    //  memset(*H, 0, size_h);
    //  memset(*C, 0, size_c);

    size_f = simd_depth * max_q_len * sizeof(float);
    if (max_q_len > *F_size) {
        if (*F_size > 0) _mm_free(*F);
        *F = (float *) _mm_malloc(size_f, 32);
        //printf("new F %x\n", *F);
        *F_size = max_q_len;
    }
    //  memset(*F, 0, size_f);

    int max_size = (max_r_len > max_q_len ? max_r_len : max_q_len);
    if (max_size > *aux_size) {
        if (*aux_size > 0) {
            free(*q_aux);
            free(*r_aux);
        }
        *q_aux = (char *) calloc(max_size * 2, sizeof(char));
        *r_aux = (char *) calloc(max_size * 2, sizeof(char));
        //printf("new q_aux %x, r_aux %x\n", *q_aux, *r_aux);
        *aux_size = max_size;
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
