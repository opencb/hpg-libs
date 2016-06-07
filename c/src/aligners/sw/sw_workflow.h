#ifndef SW_WORKFLOW_H
#define SW_WORKFLOW_H

//------------------------------------------------------------------------------------

void run_sw_workflow(char *q_filename, char *r_filename,
		     float match, float mismatch, float gap_open, float gap_extend,
		     char *matrix_filename, int batch_size, int num_threads,
		     char *out_filename);

//------------------------------------------------------------------------------------

#endif // SW_WORKFLOW_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
