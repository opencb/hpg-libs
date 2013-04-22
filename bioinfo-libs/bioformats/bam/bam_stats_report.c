/*
 * bam_stats_report.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "bam_stats_report.h"

//--------------------------------------------------------------------

void init_report_graph(report_graph_t* graph) {
    graph->x_autoscale = 1;
    graph->x_start = 1;
    graph->x_end = 100;
    graph->y_autoscale = 1;
    graph->y_start = 0;
    graph->y_end = 100;
    graph->lmargin = 10;
    graph->rmargin = 4;
    graph->tmargin = 3;
    graph->bmargin = 4;
    *graph->title = 0;
    *graph->xlabel = 0;
    *graph->ylabel = 0;
    *graph->type = 0;
    graph->x_column = 0;
    graph->num_y_columns = 1;
    graph->y_columns[0] = 1;
}

//--------------------------------------------------------------------

void generate_gnuplot_image(report_graph_t *graph, char *data_filename, char *prefix) {
    // lines specifying input data and output graph are declared and filled
    char line[2048];
    
    char gnuplot_filename[strlen(prefix) + 100];
    sprintf(gnuplot_filename, "%s.gnuplot", prefix);

    // open the file for writing gnuplot lines
    FILE* f = fopen(gnuplot_filename, "w");
    
    if (f == NULL) {
      LOG_FATAL("Impossible save file for BAM report");
    }

    sprintf(line, "set output '%s.png'\n", prefix);
    fprintf(f, line);
    fprintf(f, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
    sprintf(line, "set ylabel '%s'\n", graph->ylabel);
    fprintf(f, line);
    sprintf(line, "set xlabel '%s'\n", graph->xlabel);
    fprintf(f, line);
    fprintf(f, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
    sprintf(line, "set title '%s'\n", graph->title);
    fprintf(f, line);

    if (graph->x_autoscale == 1) {
      fprintf(f, "set autoscale x\n");
    } else {
      sprintf(line, "set xrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
      fprintf(f, line);
    }

    if (graph->y_autoscale == 1) {
      fprintf(f, "set autoscale y\n");
    } else {
      sprintf(line, "set yrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
      fprintf(f, line);
    }

    sprintf(line, "set lmargin '%i'\n", graph->lmargin);
    fprintf(f, line);
    sprintf(line, "set rmargin '%i'\n", graph->rmargin);
    fprintf(f, line);
    sprintf(line, "set tmargin '%i'\n", graph->tmargin);
    fprintf(f, line);
    sprintf(line, "set bmargin '%i'\n", graph->bmargin);
    fprintf(f, line);

    sprintf(line, "plot ");

    for (int i = 0; i < graph->num_y_columns; i++) {
      sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", 
	      line, (i == 0 ? "" : ", "), data_filename, 
	      graph->x_column, graph->y_columns[i], 
	      graph->y_titles[i], graph->type);
    }
    fprintf(f, line);
    fprintf(f, "\n");

    fclose(f);    
    
    // build the command line by calling gnuplot followed by is instruction file
    // execute command line: gnuplot filename.gnuplot
    sprintf(line, "gnuplot %s;", gnuplot_filename);
    system(line);
}

//--------------------------------------------------------------------

void report_summary(char *prefix, bam_stats_output_t *output) {

  FILE *f; 
  char path[strlen(prefix) + 100];
  sprintf(path, "%s.summary.txt", prefix);

  if ( (f = fopen(path, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  fprintf(f, "\n");
  fprintf(f, "===============================================================\n");
  fprintf(f, "=            S T A T I S T I C S     S U M M A R Y            =\n");
  fprintf(f, "===============================================================\n");
  fprintf(f, "Coverage:\n");
  fprintf(f, "\tReference length       : %lu\n", output->ref_length);
  size_t mapped_nts = output->ref_length - output->unmapped_nts;;
  fprintf(f, "\tMapped reference length: %lu (%0.2f %)\n", 
	 mapped_nts, 100.0 * mapped_nts / output->ref_length);
  fprintf(f, "\n");
  fprintf(f, "\tAverage coverage (whole reference length)      : %0.2f\n", output->depth);
  fprintf(f, "\tAverage coverage (only mapped reference length): %0.2f\n", 
	 1.0f * output->depth * output->ref_length/ mapped_nts);
  fprintf(f, "\n");
  fprintf(f, "Alignment (and read) counters:\n");
  fprintf(f, "\tNumber of alignments       : %lu\n", output->num_reads);
  fprintf(f, "\tNumber of unique alignments: %lu (%0.2f %)\n", 
	 output->num_unique_alignments, 100.0 * output->num_unique_alignments / output->num_reads);
  fprintf(f, "\n");
  fprintf(f, "\tMinimum alignment length: %lu\n", output->min_alignment_length);
  fprintf(f, "\tMaximum alignment length: %lu\n", output->max_alignment_length);
  fprintf(f, "\tAverage alignment length: %0.2f\n", 1.0f * output->num_nucleotides / output->num_reads);
  fprintf(f, "\n");
  if (output->single_end) {
    fprintf(f, "\tMode: single\n");
  } else {
    fprintf(f, "\tMode: pair\n");
  }
  fprintf(f, "\tNumber of mapped reads     : %lu (%0.2f %)\n", 
	 output->num_mapped_reads, 100.0 * output->num_mapped_reads / output->num_reads);
  fprintf(f, "\tNumber of unmapped reads   : %lu (%0.2f %)\n",
	 output->num_unmapped_reads, 100.0 * output->num_unmapped_reads / output->num_reads);
  if (!output->single_end) {
    fprintf(f, "\tNumber of mapped reads in pair #1  : %lu (%0.2f %)\n", 
	   output->num_mapped_reads_1, 100.0 * output->num_mapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #1: %lu (%0.2f %)\n", 
	   output->num_unmapped_reads_1, 100.0 * output->num_unmapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of mapped reads in pair #2  : %lu (%0.2f %)\n", 
	   output->num_mapped_reads_2, 100.0 * output->num_mapped_reads_2 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #2: %lu (%0.2f %)\n", 
	   output->num_unmapped_reads_2, 100.0 * output->num_unmapped_reads_2 / output->num_reads);
    fprintf(f, "\n");
    fprintf(f, "\tInsert size:\n");
    fprintf(f, "\t\tMinimum insert size: %lu\n", output->min_insert_size);
    fprintf(f, "\t\tMaximum insert size: %lu\n", output->max_insert_size);
    fprintf(f, "\t\tAverage insert size: %0.2f\n", 1.0f * output->insert_size_acc / output->num_reads);
  }

  fprintf(f, "\n");
  fprintf(f, "Number of errors (mismatches, indels):\n");
  for (int i = 0; i < NUM_ERRORS_STATS; i++) {
  fprintf(f, "\tNumber of reads with %i errors: %lu (%0.2f %)\n", 
	  i, output->num_errors[i], 100.0 * output->num_errors[i] / output->num_mapped_reads);
  }
  fprintf(f, "\tNumber of reads with more than %i errors: %lu (%0.2f %)\n", 
	  NUM_ERRORS_STATS, output->num_errors[NUM_ERRORS_STATS], 100.0 * output->num_errors[NUM_ERRORS_STATS] / output->num_mapped_reads);
  /*
  fprintf(f, "\tNumber of reads with 0 errors: %lu (%0.2f %)\n", 
	 output->num_0_errors, 100.0 * output->num_0_errors / output->num_mapped_reads);
  fprintf(f, "\tNumber of reads with 1 ~ 10 errors: %lu (%0.2f %)\n", 
	 output->num_1_10_errors, 100.0 * output->num_1_10_errors / output->num_mapped_reads);
  fprintf(f, "\tNumber of reads with more than 10 errors: %lu (%0.2f %)\n", 
	 output->num_11_errors, 100.0 * output->num_11_errors / output->num_mapped_reads);
  */
  fprintf(f, "\n");
  fprintf(f, "\tNumber of indels: %lu\n", 
	 output->num_indels);
  fprintf(f, "\tAverage indel size: %0.2f\n", 
	 1.0f * output->indels_acc / output->num_indels);
  fprintf(f, "\tAverage indels per read: %0.2f\n", 
	 1.0f * output->num_indels / output->num_mapped_reads);

  fprintf(f, "\n");
  fprintf(f, "Per strand:\n");
  fprintf(f, "Strand\tNum. unique alignments\tNum. mapped reads\n");
  fprintf(f, "---------------------------------------------------------------------------------------\n");
  for (int i = 0; i < 2; i++) {
    fprintf(f, "%c\t%lu (%0.2f %)\t%lu (%0.2f %)\n",
	   (i == 0 ? '+' : '-'), 
	   output->num_unique_alignments_strand[i], 100.0 * output->num_unique_alignments_strand[i] / output->num_mapped_reads,
	   output->num_mapped_reads_strand[i], 100.0 * output->num_mapped_reads_strand[i] / output->num_mapped_reads);
  }
  fprintf(f, "\n");
  fprintf(f, "Alignment quality:\n");
  fprintf(f, "\tMinimum alignment quality: %lu\n", output->min_quality);
  fprintf(f, "\tMaximum alignment quality: %lu\n", output->max_quality);
  fprintf(f, "\tAverage alignment quality: %0.2f\n", 1.0f * output->quality_acc / output->num_mapped_reads);
  fprintf(f, "\n");
  fprintf(f, "Nucleotide content:\n");
  fprintf(f, "\tNumber of A's: %lu (%0.2f %)\n", 
	 output->num_As, 100.0 * output->num_As / output->num_nucleotides);
  fprintf(f, "\tNumber of C's: %lu (%0.2f %)\n", 
	 output->num_Cs, 100.0 * output->num_Cs / output->num_nucleotides);
  fprintf(f, "\tNumber of T's: %lu (%0.2f %)\n", 
	 output->num_Ts, 100.0 * output->num_Ts / output->num_nucleotides); 
  fprintf(f, "\tNumber of G's: %lu (%0.2f %)\n", 
	 output->num_Gs, 100.0 * output->num_Gs / output->num_nucleotides);
  fprintf(f, "\tNumber of N's: %lu (%0.2f %)\n", 
	 output->num_Ns, 100.0 * output->num_Ns / output->num_nucleotides);
  fprintf(f, "\tGC percentage: %0.2f %\n", 
	 100.0 * (output->num_Cs + output->num_Gs) / output->num_nucleotides);
  fprintf(f, "\n");
  fprintf(f, "Per sequences (chromosomes) statistics:\n");
  fprintf(f, "Number of sequences: %lu\n", output->num_sequences);
  fprintf(f, "Name\tLength\tCoverage\n");
  fprintf(f, "--------------------------------------\n");
  for (int i = 0; i < output->num_sequences; i++) {
    fprintf(f, "%s\t%lu\t%0.2f\n", output->sequence_labels[i], output->sequence_lengths[i], output->depth_per_sequence[i]);
  }

  if (f != stdout) fclose(f);
}

//--------------------------------------------------------------------

void report_coverage(char *prefix, bam_stats_output_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = 100;

  // coverage histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  size_t index, len;
  for (int i = 0; i < output->num_sequences; i++) {
    len = output->sequence_lengths[i];
    for (int j = 0; j < len; j++) {
      index = output->sequence_depths_per_nt[i][j];
      if (index >= num_cols) {
	index = num_cols - 1;
      }
      hist[index]++;
    }
  }

  // coverage histogram data
  sprintf(data_filename, "%s.coverage.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i, hist[i]);
  }
  fclose(f);

  // coverage histogram image
  sprintf(img_prefix, "%s.coverage.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Coverage Histogram");
  strcpy(graph.xlabel, "Coverage");
  strcpy(graph.ylabel, "Number of genomic locations");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // coverage fraction
  double acc_hist[num_cols];
  memset(acc_hist, 0, num_cols * sizeof(double));

  size_t acc = 0;
  for (int i = num_cols - 1; i > 0; i--) {
    acc += hist[i];
    acc_hist[i] = 100.0f * acc / output->ref_length;
  }

  // coverage fraction data
  sprintf(data_filename, "%s.coverage.fraction.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%0.2f\n", i, acc_hist[i]);
  }
  fclose(f);

  // coverage histogram image
  sprintf(img_prefix, "%s.coverage.fraction", prefix);
  init_report_graph(&graph);
  
  strcpy(graph.title, "Genome Fraction Coverage");
  strcpy(graph.xlabel, "Coverage");
  strcpy(graph.ylabel, "Fraction of refence (%)");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_num_errors(char *prefix, bam_stats_output_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = NUM_ERRORS_STATS + 1;

  // num. errors histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // num. errors histogram data
  sprintf(data_filename, "%s.num.errors.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i, output->num_errors[i]);
  }
  fclose(f);

  // num. errors histogram image
  sprintf(img_prefix, "%s.num.errors.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Num. Errors Histogram");
  strcpy(graph.xlabel, "Num. errors");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_quality(char *prefix, bam_stats_output_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = QUALITY_STATS;

  // quality histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // quality histogram data
  sprintf(data_filename, "%s.quality.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i, output->quality[i]);
  }
  fclose(f);

  // quality histogram image
  sprintf(img_prefix, "%s.quality.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Quality Histogram");
  strcpy(graph.xlabel, "Quality");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_gc_content(char *prefix, bam_stats_output_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = 100;

  // GC content histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // GC content histogram data
  sprintf(data_filename, "%s.gc.content.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i + 1, output->GC_content[i]);
  }
  fclose(f);

  // GC content histogram image
  sprintf(img_prefix, "%s.gc.content.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "GC Content Distribution");
  strcpy(graph.xlabel, "GC content (%)");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_bam_stats_output(char *bam_filename, char *outdir, 
			     bam_stats_output_t *stats) {

  char prefix[strlen(bam_filename) + strlen(outdir) + 100];
  sprintf(prefix, "%s/%s", outdir, bam_filename);

  report_summary(prefix, stats);
  report_coverage(prefix, stats);
  report_num_errors(prefix, stats);
  report_quality(prefix, stats);
  report_gc_content(prefix, stats);

  //  report_insert(prefix, stats);
  //  report_nucleotides(prefix, stats);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
