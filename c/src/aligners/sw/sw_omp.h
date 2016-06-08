#ifndef SW_OMP_H
#define SW_OMP_H

//------------------------------------------------------------------------------------

void run_sw_omp(char *q_filename, char *r_filename,
		float match, float mismatch, float gap_open, float gap_extend,
		char *matrix_filename, int batch_size, int num_threads,
		char *out_filename);

//------------------------------------------------------------------------------------

#endif // SW_OMP_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
