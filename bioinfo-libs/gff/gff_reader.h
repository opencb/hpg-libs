#ifndef GFF_RAGEL_H
#define GFF_RAGEL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/linked_list.h>

#include <list.h>

#include "gff_file_structure.h"
#include "gff_file.h"
#include "gff_read.h"
#include "gff_batch.h"

enum GFF_Field { SEQUENCE, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE };

int gff_ragel_read(list_t *batches_list, size_t batch_size, gff_file_t *file);

#endif
