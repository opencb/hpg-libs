
#line 1 "ped.ragel"
#include "ped_reader.h"

static size_t lines = 0;
static size_t num_records = 0;
static size_t genotype = 0;

static ped_record_t *current_record;
static ped_batch_t *current_batch;


#line 14 "ped_reader.c"
static const int ped_start = 95;
static const int ped_first_final = 95;
static const int ped_error = 0;

static const int ped_en_main = 95;


#line 163 "ped.ragel"




int ped_ragel_read(list_t *batches_list, size_t batch_size, ped_file_t *file)
{
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_batch = ped_batch_new(batch_size);

    
#line 40 "ped_reader.c"
	{
	cs = ped_start;
	}

#line 45 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 95:
	if ( (*p) == 10 )
		goto st96;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr169;
	goto tr0;
tr0:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	goto st0;
tr3:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	goto st0;
tr7:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr11:
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr15:
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr19:
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr48:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr51:
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr55:
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr65:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr68:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr73:
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr80:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr83:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr88:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr95:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr99:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr102:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr106:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr108:
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr117:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr122:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr124:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr127:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr133:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr140:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr144:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr147:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	goto st0;
tr150:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr153:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr157:
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr195:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr221:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr224:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr227:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr252:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr273:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr278:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr281:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr299:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr322:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr328:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr346:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr351:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr356:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr374:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr377:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr380:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr385:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr388:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr391:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr398:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr403:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr408:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr416:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr421:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr426:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr431:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr436:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr441:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr446:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr451:
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr456:
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
#line 951 "ped_reader.c"
st0:
cs = 0;
	goto _out;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	if ( (*p) == 10 )
		goto st96;
	goto st0;
tr169:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 977 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr1;
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto st1;
	} else if ( (*p) >= 9 )
		goto tr1;
	goto tr0;
tr1:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 996 "ped_reader.c"
	if ( (*p) == 95 )
		goto tr4;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr4;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr4;
	} else
		goto tr4;
	goto tr3;
tr4:
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
#line 1018 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr5;
		case 95: goto st3;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr5;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st3;
		} else if ( (*p) >= 65 )
			goto st3;
	} else
		goto st3;
	goto tr3;
tr5:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 1045 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr8;
		case 95: goto tr9;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr9;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr9;
	} else
		goto tr9;
	goto tr7;
tr8:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
#line 1069 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr10;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr10;
	goto tr7;
tr10:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 1087 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr12;
		case 95: goto tr13;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr13;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr13;
	} else
		goto tr13;
	goto tr11;
tr12:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 1111 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr14;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr14;
	goto tr11;
tr14:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 1129 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr16;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr17;
	goto tr15;
tr16:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 1145 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr18;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr18;
	goto tr15;
tr18:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 1169 "ped_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr20;
	goto tr19;
tr20:
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 1183 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr21;
		case 46: goto st88;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st11;
	} else if ( (*p) >= 9 )
		goto tr21;
	goto tr19;
tr21:
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
#line 1209 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st12;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto st12;
	goto st0;
tr501:
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st97;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
#line 1240 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto st97;
	goto st0;
tr170:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st98;
tr502:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st98;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
#line 1314 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr172;
		case 65: goto tr172;
		case 67: goto tr172;
		case 71: goto tr172;
		case 78: goto tr172;
		case 84: goto tr172;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st97;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr172;
	} else
		goto tr171;
	goto tr0;
tr171:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 1352 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr27;
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto st13;
	} else if ( (*p) >= 9 )
		goto tr27;
	goto tr0;
tr27:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 1371 "ped_reader.c"
	if ( (*p) == 95 )
		goto tr29;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr3;
tr29:
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 1393 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr30;
		case 95: goto st15;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr30;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 65 )
			goto st15;
	} else
		goto st15;
	goto tr3;
tr30:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 1420 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr32;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr33;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr33;
	} else
		goto tr33;
	goto tr7;
tr32:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
#line 1444 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr34;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr34;
	goto tr7;
tr34:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 1462 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr35;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr36;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr36;
	} else
		goto tr36;
	goto tr11;
tr35:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 1486 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr37;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr37;
	goto tr11;
tr37:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 1504 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr38;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr39;
	goto tr15;
tr38:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
#line 1520 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr40;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr40;
	goto tr15;
tr40:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
#line 1544 "ped_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr41;
	goto tr19;
tr41:
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 1558 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr42;
		case 46: goto st59;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st23;
	} else if ( (*p) >= 9 )
		goto tr42;
	goto tr19;
tr42:
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
#line 1584 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto st24;
	goto st0;
tr190:
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st99;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
#line 1615 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto st99;
	goto st0;
tr173:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st100;
tr191:
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st100;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
#line 1689 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr174;
		case 65: goto tr174;
		case 67: goto tr174;
		case 71: goto tr174;
		case 78: goto tr174;
		case 84: goto tr174;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st99;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr174;
	} else
		goto tr171;
	goto tr0;
tr174:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
	goto st101;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
#line 1727 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr176;
		case 32: goto tr175;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr175;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr0;
tr175:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st102;
tr197:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st102;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
#line 1775 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 66: goto tr29;
		case 71: goto tr178;
		case 78: goto tr178;
		case 84: goto tr178;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr178;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr178;
	goto tr3;
tr178:
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st103;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
#line 1811 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr180;
		case 32: goto tr179;
		case 46: goto st99;
		case 66: goto st15;
		case 71: goto st103;
		case 78: goto st103;
		case 84: goto st103;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st103;
		} else if ( (*p) >= 9 )
			goto tr179;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st103;
	goto tr3;
tr179:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st104;
tr229:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st104;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
#line 1862 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr182;
		case 66: goto tr33;
		case 71: goto tr183;
		case 78: goto tr183;
		case 84: goto tr183;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr183;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr183;
	goto tr7;
tr182:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st105;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
#line 1898 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr185;
		case 32: goto tr184;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr184;
	goto tr7;
tr184:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st106;
tr287:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
#line 1944 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr138;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr94;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr11;
tr138:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
#line 1980 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr187;
		case 32: goto tr186;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr186;
	goto tr11;
tr186:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st108;
tr245:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st108;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
#line 2026 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr77;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr112;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr15;
tr77:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st109;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
#line 2053 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto tr188;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr188;
	goto tr15;
tr188:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st110;
tr216:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st110;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
#line 2111 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr60;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr19;
tr60:
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st111;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
#line 2138 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr191;
		case 32: goto tr190;
		case 46: goto st112;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st111;
	} else if ( (*p) >= 9 )
		goto tr190;
	goto tr19;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st113;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr19;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 10: goto tr191;
		case 32: goto tr190;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st113;
	} else if ( (*p) >= 9 )
		goto tr190;
	goto tr19;
tr189:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st114;
tr217:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st114;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
#line 2276 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr174;
		case 65: goto tr174;
		case 67: goto tr174;
		case 71: goto tr174;
		case 78: goto tr174;
		case 84: goto tr174;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st99;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr196;
	} else
		goto tr171;
	goto tr195;
tr196:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st115;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
#line 2318 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr198;
		case 32: goto tr197;
		case 46: goto st236;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr197;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st115;
	} else
		goto st13;
	goto tr195;
tr176:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st116;
tr198:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st116;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
#line 2406 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr174;
		case 66: goto tr202;
		case 71: goto tr201;
		case 78: goto tr201;
		case 84: goto tr201;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr201;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr201;
	goto tr147;
tr201:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st117;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
#line 2463 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr204;
		case 32: goto tr203;
		case 46: goto st101;
		case 66: goto st79;
		case 71: goto st117;
		case 78: goto st117;
		case 84: goto st117;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr203;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st117;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st117;
	goto tr147;
tr203:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st118;
tr254:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st118;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
#line 2534 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr182;
		case 66: goto tr105;
		case 71: goto tr206;
		case 78: goto tr206;
		case 84: goto tr206;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr206;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr206;
	goto tr122;
tr206:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st119;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
#line 2574 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr208;
		case 32: goto tr207;
		case 46: goto st99;
		case 66: goto st65;
		case 71: goto st119;
		case 78: goto st119;
		case 84: goto st119;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st119;
		} else if ( (*p) >= 9 )
			goto tr207;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st119;
	goto tr122;
tr207:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st120;
tr330:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st120;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
#line 2637 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr210;
		case 66: goto tr71;
		case 71: goto tr211;
		case 78: goto tr211;
		case 84: goto tr211;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr211;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr211;
	goto tr48;
tr210:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st121;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
#line 2677 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr213;
		case 32: goto tr212;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr212;
	goto tr48;
tr212:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st122;
tr315:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st122;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
#line 2735 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr92;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr121;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr51;
tr92:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st123;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
#line 2775 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr215;
		case 32: goto tr214;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr214;
	goto tr51;
tr214:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st124;
tr270:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st124;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
#line 2845 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr77;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr78;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr55;
tr78:
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st125;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
#line 2876 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr217;
		case 32: goto tr216;
		case 46: goto st126;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st125;
	} else if ( (*p) >= 9 )
		goto tr216;
	goto tr55;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st127;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr55;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	switch( (*p) ) {
		case 10: goto tr217;
		case 32: goto tr216;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st127;
	} else if ( (*p) >= 9 )
		goto tr216;
	goto tr55;
tr215:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st128;
tr271:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st128;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
#line 3026 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr222;
		case 65: goto tr174;
		case 67: goto tr174;
		case 71: goto tr174;
		case 78: goto tr174;
		case 84: goto tr174;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st99;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr223;
	} else
		goto tr171;
	goto tr221;
tr222:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st129;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
#line 3068 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr226;
		case 32: goto tr225;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr225;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr224;
tr225:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st130;
tr462:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st130;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
#line 3140 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 66: goto tr29;
		case 71: goto tr178;
		case 78: goto tr178;
		case 84: goto tr178;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr228;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr178;
	goto tr227;
tr228:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st131;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
#line 3180 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr230;
		case 32: goto tr229;
		case 46: goto st112;
		case 66: goto st15;
		case 71: goto st103;
		case 78: goto st103;
		case 84: goto st103;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st131;
		} else if ( (*p) >= 9 )
			goto tr229;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st103;
	goto tr227;
tr180:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st132;
tr230:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st132;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
#line 3271 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr232;
		case 66: goto tr234;
		case 71: goto tr233;
		case 78: goto tr233;
		case 84: goto tr233;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr233;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr233;
	goto tr144;
tr232:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st133;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
#line 3328 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr236;
		case 32: goto tr235;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr235;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr144;
tr235:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st134;
tr301:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st134;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
#line 3388 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr138;
		case 66: goto tr129;
		case 71: goto tr237;
		case 78: goto tr237;
		case 84: goto tr237;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr237;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr237;
	goto tr140;
tr237:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st135;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
#line 3428 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr239;
		case 32: goto tr238;
		case 46: goto st99;
		case 66: goto st74;
		case 71: goto st135;
		case 78: goto st135;
		case 84: goto st135;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st135;
		} else if ( (*p) >= 9 )
			goto tr238;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st135;
	goto tr140;
tr238:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st136;
tr393:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st136;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
#line 3491 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr241;
		case 66: goto tr33;
		case 71: goto tr183;
		case 78: goto tr183;
		case 84: goto tr183;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr242;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr183;
	goto tr106;
tr241:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st137;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
#line 3531 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr244;
		case 32: goto tr243;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr243;
	goto tr106;
tr243:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st138;
tr343:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st138;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
#line 3601 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr138;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr139;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr108;
tr139:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st139;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
#line 3641 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr246;
		case 32: goto tr245;
		case 46: goto st112;
		case 66: goto st25;
		case 71: goto st152;
		case 78: goto st152;
		case 84: goto st152;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st139;
		} else if ( (*p) >= 9 )
			goto tr245;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st152;
	goto tr108;
tr187:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st140;
tr246:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st140;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
#line 3736 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr222;
		case 65: goto tr174;
		case 67: goto tr174;
		case 71: goto tr174;
		case 78: goto tr174;
		case 84: goto tr174;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st99;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr249;
	} else
		goto tr171;
	goto tr224;
tr249:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st141;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
#line 3778 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr226;
		case 32: goto tr225;
		case 46: goto st234;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr225;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st141;
	} else
		goto st13;
	goto tr224;
tr226:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st142;
tr463:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st142;
st142:
	if ( ++p == pe )
		goto _test_eof142;
case 142:
#line 3890 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr174;
		case 66: goto tr202;
		case 71: goto tr201;
		case 78: goto tr201;
		case 84: goto tr201;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr253;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr201;
	goto tr252;
tr253:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st143;
st143:
	if ( ++p == pe )
		goto _test_eof143;
case 143:
#line 3951 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr255;
		case 32: goto tr254;
		case 46: goto st236;
		case 66: goto st79;
		case 71: goto st117;
		case 78: goto st117;
		case 84: goto st117;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr254;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st143;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st117;
	goto tr252;
tr204:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st144;
tr255:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st144;
st144:
	if ( ++p == pe )
		goto _test_eof144;
case 144:
#line 4062 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr232;
		case 66: goto tr258;
		case 71: goto tr257;
		case 78: goto tr257;
		case 84: goto tr257;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr257;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr257;
	goto tr95;
tr257:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st145;
st145:
	if ( ++p == pe )
		goto _test_eof145;
case 145:
#line 4123 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr260;
		case 32: goto tr259;
		case 46: goto st101;
		case 66: goto st50;
		case 71: goto st145;
		case 78: goto st145;
		case 84: goto st145;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr259;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st145;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st145;
	goto tr95;
tr259:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st146;
tr358:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st146;
st146:
	if ( ++p == pe )
		goto _test_eof146;
case 146:
#line 4206 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr210;
		case 66: goto tr85;
		case 71: goto tr262;
		case 78: goto tr262;
		case 84: goto tr262;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr262;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr262;
	goto tr65;
tr262:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st147;
st147:
	if ( ++p == pe )
		goto _test_eof147;
case 147:
#line 4250 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr264;
		case 32: goto tr263;
		case 46: goto st99;
		case 66: goto st37;
		case 71: goto st147;
		case 78: goto st147;
		case 84: goto st147;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st147;
		} else if ( (*p) >= 9 )
			goto tr263;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st147;
	goto tr65;
tr263:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st148;
tr410:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st148;
st148:
	if ( ++p == pe )
		goto _test_eof148;
case 148:
#line 4325 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr266;
		case 66: goto tr71;
		case 71: goto tr211;
		case 78: goto tr211;
		case 84: goto tr211;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr267;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr211;
	goto tr68;
tr266:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st149;
st149:
	if ( ++p == pe )
		goto _test_eof149;
case 149:
#line 4369 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr269;
		case 32: goto tr268;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st99;
	} else if ( (*p) >= 9 )
		goto tr268;
	goto tr68;
tr268:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st150;
tr371:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st150;
st150:
	if ( ++p == pe )
		goto _test_eof150;
case 150:
#line 4451 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr92;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr93;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr73;
tr93:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st151;
st151:
	if ( ++p == pe )
		goto _test_eof151;
case 151:
#line 4495 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr271;
		case 32: goto tr270;
		case 46: goto st126;
		case 66: goto st25;
		case 71: goto st152;
		case 78: goto st152;
		case 84: goto st152;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st151;
		} else if ( (*p) >= 9 )
			goto tr270;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st152;
	goto tr73;
tr94:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st152;
st152:
	if ( ++p == pe )
		goto _test_eof152;
case 152:
#line 4531 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr187;
		case 32: goto tr186;
		case 46: goto st99;
		case 66: goto st25;
		case 71: goto st152;
		case 78: goto st152;
		case 84: goto st152;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st152;
		} else if ( (*p) >= 9 )
			goto tr186;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st152;
	goto tr11;
tr36:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 4567 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr37;
		case 95: goto st25;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr37;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 65 )
			goto st25;
	} else
		goto st25;
	goto tr11;
tr269:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st153;
tr372:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st153;
st153:
	if ( ++p == pe )
		goto _test_eof153;
case 153:
#line 4689 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr274;
		case 66: goto tr277;
		case 71: goto tr276;
		case 78: goto tr276;
		case 84: goto tr276;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr275;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr276;
	goto tr273;
tr274:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st154;
st154:
	if ( ++p == pe )
		goto _test_eof154;
case 154:
#line 4750 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr280;
		case 32: goto tr279;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr279;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr278;
tr279:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st155;
tr475:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st155;
st155:
	if ( ++p == pe )
		goto _test_eof155;
case 155:
#line 4834 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr77;
		case 66: goto tr29;
		case 71: goto tr178;
		case 78: goto tr178;
		case 84: goto tr178;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr282;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr178;
	goto tr281;
tr282:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st156;
st156:
	if ( ++p == pe )
		goto _test_eof156;
case 156:
#line 4878 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr284;
		case 32: goto tr283;
		case 46: goto st126;
		case 66: goto st15;
		case 71: goto st103;
		case 78: goto st103;
		case 84: goto st103;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st156;
		} else if ( (*p) >= 9 )
			goto tr283;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st103;
	goto tr281;
tr295:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st157;
tr283:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st157;
st157:
	if ( ++p == pe )
		goto _test_eof157;
case 157:
#line 4953 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr182;
		case 66: goto tr33;
		case 71: goto tr183;
		case 78: goto tr183;
		case 84: goto tr183;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr286;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr183;
	goto tr157;
tr286:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st158;
st158:
	if ( ++p == pe )
		goto _test_eof158;
case 158:
#line 4993 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr288;
		case 32: goto tr287;
		case 46: goto st112;
		case 66: goto st73;
		case 71: goto st247;
		case 78: goto st247;
		case 84: goto st247;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st158;
		} else if ( (*p) >= 9 )
			goto tr287;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st247;
	goto tr157;
tr185:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st159;
tr288:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st159;
st159:
	if ( ++p == pe )
		goto _test_eof159;
case 159:
#line 5088 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr291;
		case 66: goto tr277;
		case 71: goto tr276;
		case 78: goto tr276;
		case 84: goto tr276;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr276;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr276;
	goto tr150;
tr291:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st160;
st160:
	if ( ++p == pe )
		goto _test_eof160;
case 160:
#line 5145 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr293;
		case 32: goto tr292;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr292;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr150;
tr292:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st161;
tr348:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st161;
st161:
	if ( ++p == pe )
		goto _test_eof161;
case 161:
#line 5205 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr77;
		case 66: goto tr29;
		case 71: goto tr178;
		case 78: goto tr178;
		case 84: goto tr178;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr294;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr178;
	goto tr153;
tr294:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st162;
st162:
	if ( ++p == pe )
		goto _test_eof162;
case 162:
#line 5245 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr296;
		case 32: goto tr295;
		case 46: goto st239;
		case 66: goto st15;
		case 71: goto st103;
		case 78: goto st103;
		case 84: goto st103;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st162;
		} else if ( (*p) >= 9 )
			goto tr295;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st103;
	goto tr153;
tr284:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st163;
tr296:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st163;
st163:
	if ( ++p == pe )
		goto _test_eof163;
case 163:
#line 5360 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr232;
		case 66: goto tr234;
		case 71: goto tr233;
		case 78: goto tr233;
		case 84: goto tr233;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr300;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr233;
	goto tr299;
tr300:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st164;
st164:
	if ( ++p == pe )
		goto _test_eof164;
case 164:
#line 5421 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr302;
		case 32: goto tr301;
		case 46: goto st236;
		case 66: goto st77;
		case 71: goto st243;
		case 78: goto st243;
		case 84: goto st243;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr301;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st164;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st243;
	goto tr299;
tr236:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st165;
tr302:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st165;
st165:
	if ( ++p == pe )
		goto _test_eof165;
case 165:
#line 5536 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr291;
		case 66: goto tr306;
		case 71: goto tr305;
		case 78: goto tr305;
		case 84: goto tr305;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr305;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr305;
	goto tr99;
tr305:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st166;
st166:
	if ( ++p == pe )
		goto _test_eof166;
case 166:
#line 5597 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr308;
		case 32: goto tr307;
		case 46: goto st101;
		case 66: goto st53;
		case 71: goto st166;
		case 78: goto st166;
		case 84: goto st166;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr307;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st166;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st166;
	goto tr99;
tr307:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st167;
tr400:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st167;
st167:
	if ( ++p == pe )
		goto _test_eof167;
case 167:
#line 5680 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr241;
		case 66: goto tr105;
		case 71: goto tr206;
		case 78: goto tr206;
		case 84: goto tr206;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr310;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr206;
	goto tr102;
tr310:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st168;
st168:
	if ( ++p == pe )
		goto _test_eof168;
case 168:
#line 5724 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr312;
		case 32: goto tr311;
		case 46: goto st239;
		case 66: goto st65;
		case 71: goto st119;
		case 78: goto st119;
		case 84: goto st119;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st168;
		} else if ( (*p) >= 9 )
			goto tr311;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st119;
	goto tr102;
tr311:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st169;
tr423:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st169;
st169:
	if ( ++p == pe )
		goto _test_eof169;
case 169:
#line 5811 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr210;
		case 66: goto tr71;
		case 71: goto tr211;
		case 78: goto tr211;
		case 84: goto tr211;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr314;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr211;
	goto tr117;
tr314:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st170;
st170:
	if ( ++p == pe )
		goto _test_eof170;
case 170:
#line 5855 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr316;
		case 32: goto tr315;
		case 46: goto st112;
		case 66: goto st26;
		case 71: goto st195;
		case 78: goto st195;
		case 84: goto st195;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st170;
		} else if ( (*p) >= 9 )
			goto tr315;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st195;
	goto tr117;
tr213:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st171;
tr316:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st171;
st171:
	if ( ++p == pe )
		goto _test_eof171;
case 171:
#line 5962 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr274;
		case 66: goto tr277;
		case 71: goto tr276;
		case 78: goto tr276;
		case 84: goto tr276;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr319;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr276;
	goto tr278;
tr319:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st172;
st172:
	if ( ++p == pe )
		goto _test_eof172;
case 172:
#line 6023 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr280;
		case 32: goto tr279;
		case 46: goto st234;
		case 66: goto st81;
		case 71: goto st246;
		case 78: goto st246;
		case 84: goto st246;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr279;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st172;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st246;
	goto tr278;
tr280:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st173;
tr476:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st173;
st173:
	if ( ++p == pe )
		goto _test_eof173;
case 173:
#line 6162 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr222;
		case 66: goto tr202;
		case 71: goto tr201;
		case 78: goto tr201;
		case 84: goto tr201;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr323;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr201;
	goto tr322;
tr323:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st174;
st174:
	if ( ++p == pe )
		goto _test_eof174;
case 174:
#line 6227 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr325;
		case 32: goto tr324;
		case 46: goto st232;
		case 66: goto st79;
		case 71: goto st117;
		case 78: goto st117;
		case 84: goto st117;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr324;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st174;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st117;
	goto tr322;
tr324:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st175;
tr353:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st175;
st175:
	if ( ++p == pe )
		goto _test_eof175;
case 175:
#line 6322 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr182;
		case 66: goto tr105;
		case 71: goto tr206;
		case 78: goto tr206;
		case 84: goto tr206;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr329;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr206;
	goto tr328;
tr329:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st176;
st176:
	if ( ++p == pe )
		goto _test_eof176;
case 176:
#line 6366 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr331;
		case 32: goto tr330;
		case 46: goto st112;
		case 66: goto st65;
		case 71: goto st119;
		case 78: goto st119;
		case 84: goto st119;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st176;
		} else if ( (*p) >= 9 )
			goto tr330;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st119;
	goto tr328;
tr208:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st177;
tr331:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st177;
st177:
	if ( ++p == pe )
		goto _test_eof177;
case 177:
#line 6469 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr333;
		case 66: goto tr335;
		case 71: goto tr334;
		case 78: goto tr334;
		case 84: goto tr334;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr334;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr334;
	goto tr124;
tr333:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st178;
st178:
	if ( ++p == pe )
		goto _test_eof178;
case 178:
#line 6530 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr337;
		case 32: goto tr336;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr336;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr124;
tr336:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st179;
tr428:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st179;
st179:
	if ( ++p == pe )
		goto _test_eof179;
case 179:
#line 6602 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr92;
		case 66: goto tr129;
		case 71: goto tr237;
		case 78: goto tr237;
		case 84: goto tr237;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr338;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr237;
	goto tr127;
tr338:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st180;
st180:
	if ( ++p == pe )
		goto _test_eof180;
case 180:
#line 6646 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr340;
		case 32: goto tr339;
		case 46: goto st239;
		case 66: goto st74;
		case 71: goto st135;
		case 78: goto st135;
		case 84: goto st135;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st180;
		} else if ( (*p) >= 9 )
			goto tr339;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st135;
	goto tr127;
tr339:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st181;
tr382:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st181;
st181:
	if ( ++p == pe )
		goto _test_eof181;
case 181:
#line 6733 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr241;
		case 66: goto tr33;
		case 71: goto tr183;
		case 78: goto tr183;
		case 84: goto tr183;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr342;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr183;
	goto tr133;
tr342:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st182;
st182:
	if ( ++p == pe )
		goto _test_eof182;
case 182:
#line 6777 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr344;
		case 32: goto tr343;
		case 46: goto st126;
		case 66: goto st73;
		case 71: goto st247;
		case 78: goto st247;
		case 84: goto st247;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st182;
		} else if ( (*p) >= 9 )
			goto tr343;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st247;
	goto tr133;
tr244:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st183;
tr344:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st183;
st183:
	if ( ++p == pe )
		goto _test_eof183;
case 183:
#line 6896 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr291;
		case 66: goto tr277;
		case 71: goto tr276;
		case 78: goto tr276;
		case 84: goto tr276;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr347;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr276;
	goto tr346;
tr347:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st184;
st184:
	if ( ++p == pe )
		goto _test_eof184;
case 184:
#line 6957 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr349;
		case 32: goto tr348;
		case 46: goto st236;
		case 66: goto st81;
		case 71: goto st246;
		case 78: goto st246;
		case 84: goto st246;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr348;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st184;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st246;
	goto tr346;
tr293:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st185;
tr349:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st185;
st185:
	if ( ++p == pe )
		goto _test_eof185;
case 185:
#line 7072 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr222;
		case 66: goto tr202;
		case 71: goto tr201;
		case 78: goto tr201;
		case 84: goto tr201;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr352;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr201;
	goto tr351;
tr352:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st186;
st186:
	if ( ++p == pe )
		goto _test_eof186;
case 186:
#line 7133 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr354;
		case 32: goto tr353;
		case 46: goto st234;
		case 66: goto st79;
		case 71: goto st117;
		case 78: goto st117;
		case 84: goto st117;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr353;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st186;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st117;
	goto tr351;
tr325:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st187;
tr354:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st187;
st187:
	if ( ++p == pe )
		goto _test_eof187;
case 187:
#line 7268 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr232;
		case 66: goto tr258;
		case 71: goto tr257;
		case 78: goto tr257;
		case 84: goto tr257;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr357;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr257;
	goto tr356;
tr357:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st188;
st188:
	if ( ++p == pe )
		goto _test_eof188;
case 188:
#line 7333 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr359;
		case 32: goto tr358;
		case 46: goto st236;
		case 66: goto st50;
		case 71: goto st145;
		case 78: goto st145;
		case 84: goto st145;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr358;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st188;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st145;
	goto tr356;
tr260:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st189;
tr359:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st189;
st189:
	if ( ++p == pe )
		goto _test_eof189;
case 189:
#line 7456 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr333;
		case 66: goto tr362;
		case 71: goto tr361;
		case 78: goto tr361;
		case 84: goto tr361;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr361;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr361;
	goto tr80;
tr361:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st190;
st190:
	if ( ++p == pe )
		goto _test_eof190;
case 190:
#line 7521 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr364;
		case 32: goto tr363;
		case 46: goto st101;
		case 66: goto st44;
		case 71: goto st190;
		case 78: goto st190;
		case 84: goto st190;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr363;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st190;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st190;
	goto tr80;
tr363:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st191;
tr443:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st191;
st191:
	if ( ++p == pe )
		goto _test_eof191;
case 191:
#line 7616 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr266;
		case 66: goto tr85;
		case 71: goto tr262;
		case 78: goto tr262;
		case 84: goto tr262;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr366;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr262;
	goto tr83;
tr366:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st192;
st192:
	if ( ++p == pe )
		goto _test_eof192;
case 192:
#line 7664 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr368;
		case 32: goto tr367;
		case 46: goto st239;
		case 66: goto st37;
		case 71: goto st147;
		case 78: goto st147;
		case 84: goto st147;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st192;
		} else if ( (*p) >= 9 )
			goto tr367;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st147;
	goto tr83;
tr367:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st193;
tr453:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st193;
st193:
	if ( ++p == pe )
		goto _test_eof193;
case 193:
#line 7763 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr266;
		case 66: goto tr71;
		case 71: goto tr211;
		case 78: goto tr211;
		case 84: goto tr211;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr370;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr211;
	goto tr88;
tr370:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st194;
st194:
	if ( ++p == pe )
		goto _test_eof194;
case 194:
#line 7811 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr372;
		case 32: goto tr371;
		case 46: goto st126;
		case 66: goto st26;
		case 71: goto st195;
		case 78: goto st195;
		case 84: goto st195;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st194;
		} else if ( (*p) >= 9 )
			goto tr371;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st195;
	goto tr88;
tr211:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st195;
st195:
	if ( ++p == pe )
		goto _test_eof195;
case 195:
#line 7851 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr213;
		case 32: goto tr212;
		case 46: goto st99;
		case 66: goto st26;
		case 71: goto st195;
		case 78: goto st195;
		case 84: goto st195;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st195;
		} else if ( (*p) >= 9 )
			goto tr212;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st195;
	goto tr48;
tr71:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 7891 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr49;
		case 95: goto st26;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr49;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 65 )
			goto st26;
	} else
		goto st26;
	goto tr48;
tr49:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 7926 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr52;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr53;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr36;
	} else
		goto tr36;
	goto tr51;
tr52:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st28;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
#line 7954 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr54;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr54;
	goto tr51;
tr54:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st29;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
#line 7984 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr38;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr56;
	goto tr55;
tr56:
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st30;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
#line 8004 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr57;
		case 46: goto st32;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st30;
	} else if ( (*p) >= 9 )
		goto tr57;
	goto tr55;
tr57:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
#line 8042 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr60;
	} else if ( (*p) >= 9 )
		goto st24;
	goto tr19;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st33;
	goto tr55;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 32 )
		goto tr57;
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st33;
	} else if ( (*p) >= 9 )
		goto tr57;
	goto tr55;
tr53:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st34;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
#line 8091 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr54;
		case 46: goto st35;
		case 95: goto st25;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr54;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 65 )
			goto st25;
	} else
		goto st34;
	goto tr51;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st36;
	goto tr15;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	if ( (*p) == 32 )
		goto tr40;
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st36;
	} else if ( (*p) >= 9 )
		goto tr40;
	goto tr15;
tr368:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st196;
tr454:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st196;
st196:
	if ( ++p == pe )
		goto _test_eof196;
case 196:
#line 8241 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr375;
		case 66: goto tr335;
		case 71: goto tr334;
		case 78: goto tr334;
		case 84: goto tr334;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr376;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr334;
	goto tr374;
tr375:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st197;
st197:
	if ( ++p == pe )
		goto _test_eof197;
case 197:
#line 8306 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr379;
		case 32: goto tr378;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr378;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr377;
tr378:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st198;
tr472:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st198;
st198:
	if ( ++p == pe )
		goto _test_eof198;
case 198:
#line 8402 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr92;
		case 66: goto tr129;
		case 71: goto tr237;
		case 78: goto tr237;
		case 84: goto tr237;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr381;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr237;
	goto tr380;
tr381:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st199;
st199:
	if ( ++p == pe )
		goto _test_eof199;
case 199:
#line 8450 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr383;
		case 32: goto tr382;
		case 46: goto st126;
		case 66: goto st74;
		case 71: goto st135;
		case 78: goto st135;
		case 84: goto st135;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st199;
		} else if ( (*p) >= 9 )
			goto tr382;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st135;
	goto tr380;
tr340:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st200;
tr383:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st200;
st200:
	if ( ++p == pe )
		goto _test_eof200;
case 200:
#line 8577 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr386;
		case 66: goto tr234;
		case 71: goto tr233;
		case 78: goto tr233;
		case 84: goto tr233;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr387;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr233;
	goto tr385;
tr386:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st201;
st201:
	if ( ++p == pe )
		goto _test_eof201;
case 201:
#line 8638 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr390;
		case 32: goto tr389;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr389;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st101;
	} else
		goto st13;
	goto tr388;
tr389:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st202;
tr469:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st202;
st202:
	if ( ++p == pe )
		goto _test_eof202;
case 202:
#line 8722 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr138;
		case 66: goto tr129;
		case 71: goto tr237;
		case 78: goto tr237;
		case 84: goto tr237;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr392;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr237;
	goto tr391;
tr392:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st203;
st203:
	if ( ++p == pe )
		goto _test_eof203;
case 203:
#line 8766 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr394;
		case 32: goto tr393;
		case 46: goto st112;
		case 66: goto st74;
		case 71: goto st135;
		case 78: goto st135;
		case 84: goto st135;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st203;
		} else if ( (*p) >= 9 )
			goto tr393;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st135;
	goto tr391;
tr239:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st204;
tr394:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st204;
st204:
	if ( ++p == pe )
		goto _test_eof204;
case 204:
#line 8869 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr386;
		case 66: goto tr234;
		case 71: goto tr233;
		case 78: goto tr233;
		case 84: goto tr233;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr396;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr233;
	goto tr388;
tr396:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st205;
st205:
	if ( ++p == pe )
		goto _test_eof205;
case 205:
#line 8930 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr390;
		case 32: goto tr389;
		case 46: goto st234;
		case 66: goto st77;
		case 71: goto st243;
		case 78: goto st243;
		case 84: goto st243;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr389;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st205;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st243;
	goto tr388;
tr390:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st206;
tr470:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st206;
st206:
	if ( ++p == pe )
		goto _test_eof206;
case 206:
#line 9069 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr291;
		case 66: goto tr306;
		case 71: goto tr305;
		case 78: goto tr305;
		case 84: goto tr305;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr399;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr305;
	goto tr398;
tr399:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st207;
st207:
	if ( ++p == pe )
		goto _test_eof207;
case 207:
#line 9134 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr401;
		case 32: goto tr400;
		case 46: goto st236;
		case 66: goto st53;
		case 71: goto st166;
		case 78: goto st166;
		case 84: goto st166;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr400;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st207;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st166;
	goto tr398;
tr308:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st208;
tr401:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st208;
st208:
	if ( ++p == pe )
		goto _test_eof208;
case 208:
#line 9257 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr386;
		case 66: goto tr258;
		case 71: goto tr257;
		case 78: goto tr257;
		case 84: goto tr257;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr404;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr257;
	goto tr403;
tr404:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st209;
st209:
	if ( ++p == pe )
		goto _test_eof209;
case 209:
#line 9322 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr406;
		case 32: goto tr405;
		case 46: goto st234;
		case 66: goto st50;
		case 71: goto st145;
		case 78: goto st145;
		case 84: goto st145;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr405;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st209;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st145;
	goto tr403;
tr405:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st210;
tr438:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st210;
st210:
	if ( ++p == pe )
		goto _test_eof210;
case 210:
#line 9429 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr210;
		case 66: goto tr85;
		case 71: goto tr262;
		case 78: goto tr262;
		case 84: goto tr262;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr409;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr262;
	goto tr408;
tr409:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st211;
st211:
	if ( ++p == pe )
		goto _test_eof211;
case 211:
#line 9477 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr411;
		case 32: goto tr410;
		case 46: goto st112;
		case 66: goto st37;
		case 71: goto st147;
		case 78: goto st147;
		case 84: goto st147;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st211;
		} else if ( (*p) >= 9 )
			goto tr410;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st147;
	goto tr408;
tr264:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st212;
tr411:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st212;
st212:
	if ( ++p == pe )
		goto _test_eof212;
case 212:
#line 9592 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr375;
		case 66: goto tr335;
		case 71: goto tr334;
		case 78: goto tr334;
		case 84: goto tr334;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr413;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr334;
	goto tr377;
tr413:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st213;
st213:
	if ( ++p == pe )
		goto _test_eof213;
case 213:
#line 9657 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr379;
		case 32: goto tr378;
		case 46: goto st234;
		case 66: goto st67;
		case 71: goto st242;
		case 78: goto st242;
		case 84: goto st242;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr378;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st213;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st242;
	goto tr377;
tr379:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st214;
tr473:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st214;
st214:
	if ( ++p == pe )
		goto _test_eof214;
case 214:
#line 9808 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr274;
		case 66: goto tr306;
		case 71: goto tr305;
		case 78: goto tr305;
		case 84: goto tr305;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr417;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr305;
	goto tr416;
tr417:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st215;
st215:
	if ( ++p == pe )
		goto _test_eof215;
case 215:
#line 9877 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr419;
		case 32: goto tr418;
		case 46: goto st232;
		case 66: goto st53;
		case 71: goto st166;
		case 78: goto st166;
		case 84: goto st166;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr418;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st215;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st166;
	goto tr416;
tr418:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st216;
tr433:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st216;
st216:
	if ( ++p == pe )
		goto _test_eof216;
case 216:
#line 9984 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr241;
		case 66: goto tr105;
		case 71: goto tr206;
		case 78: goto tr206;
		case 84: goto tr206;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr422;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr206;
	goto tr421;
tr422:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st217;
st217:
	if ( ++p == pe )
		goto _test_eof217;
case 217:
#line 10032 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr424;
		case 32: goto tr423;
		case 46: goto st126;
		case 66: goto st65;
		case 71: goto st119;
		case 78: goto st119;
		case 84: goto st119;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st217;
		} else if ( (*p) >= 9 )
			goto tr423;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st119;
	goto tr421;
tr312:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st218;
tr424:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st218;
st218:
	if ( ++p == pe )
		goto _test_eof218;
case 218:
#line 10159 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr333;
		case 66: goto tr335;
		case 71: goto tr334;
		case 78: goto tr334;
		case 84: goto tr334;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr427;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr334;
	goto tr426;
tr427:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st219;
st219:
	if ( ++p == pe )
		goto _test_eof219;
case 219:
#line 10224 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr429;
		case 32: goto tr428;
		case 46: goto st236;
		case 66: goto st67;
		case 71: goto st242;
		case 78: goto st242;
		case 84: goto st242;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr428;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st219;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st242;
	goto tr426;
tr337:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st220;
tr429:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st220;
st220:
	if ( ++p == pe )
		goto _test_eof220;
case 220:
#line 10351 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr274;
		case 66: goto tr306;
		case 71: goto tr305;
		case 78: goto tr305;
		case 84: goto tr305;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr432;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr305;
	goto tr431;
tr432:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st221;
st221:
	if ( ++p == pe )
		goto _test_eof221;
case 221:
#line 10416 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr434;
		case 32: goto tr433;
		case 46: goto st234;
		case 66: goto st53;
		case 71: goto st166;
		case 78: goto st166;
		case 84: goto st166;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr433;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st221;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st166;
	goto tr431;
tr419:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st222;
tr434:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st222;
st222:
	if ( ++p == pe )
		goto _test_eof222;
case 222:
#line 10563 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr386;
		case 66: goto tr258;
		case 71: goto tr257;
		case 78: goto tr257;
		case 84: goto tr257;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr437;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr257;
	goto tr436;
tr437:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st223;
st223:
	if ( ++p == pe )
		goto _test_eof223;
case 223:
#line 10632 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr439;
		case 32: goto tr438;
		case 46: goto st232;
		case 66: goto st50;
		case 71: goto st145;
		case 78: goto st145;
		case 84: goto st145;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr438;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st223;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st145;
	goto tr436;
tr406:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st224;
tr439:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st224;
st224:
	if ( ++p == pe )
		goto _test_eof224;
case 224:
#line 10779 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr333;
		case 66: goto tr362;
		case 71: goto tr361;
		case 78: goto tr361;
		case 84: goto tr361;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr442;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr361;
	goto tr441;
tr442:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st225;
st225:
	if ( ++p == pe )
		goto _test_eof225;
case 225:
#line 10848 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr444;
		case 32: goto tr443;
		case 46: goto st236;
		case 66: goto st44;
		case 71: goto st190;
		case 78: goto st190;
		case 84: goto st190;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr443;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st225;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st190;
	goto tr441;
tr364:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st226;
tr444:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st226;
st226:
	if ( ++p == pe )
		goto _test_eof226;
case 226:
#line 10983 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr375;
		case 66: goto tr362;
		case 71: goto tr361;
		case 78: goto tr361;
		case 84: goto tr361;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr447;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr361;
	goto tr446;
tr447:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st227;
st227:
	if ( ++p == pe )
		goto _test_eof227;
case 227:
#line 11052 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr449;
		case 32: goto tr448;
		case 46: goto st234;
		case 66: goto st44;
		case 71: goto st190;
		case 78: goto st190;
		case 84: goto st190;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr448;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st227;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st190;
	goto tr446;
tr448:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st228;
tr458:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st228;
st228:
	if ( ++p == pe )
		goto _test_eof228;
case 228:
#line 11171 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr266;
		case 66: goto tr85;
		case 71: goto tr262;
		case 78: goto tr262;
		case 84: goto tr262;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr452;
		} else if ( (*p) >= 9 )
			goto st99;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr262;
	goto tr451;
tr452:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st229;
st229:
	if ( ++p == pe )
		goto _test_eof229;
case 229:
#line 11223 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr454;
		case 32: goto tr453;
		case 46: goto st126;
		case 66: goto st37;
		case 71: goto st147;
		case 78: goto st147;
		case 84: goto st147;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st229;
		} else if ( (*p) >= 9 )
			goto tr453;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st147;
	goto tr451;
tr85:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 11267 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr66;
		case 95: goto st37;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr66;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 65 )
			goto st37;
	} else
		goto st37;
	goto tr65;
tr66:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 11306 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr69;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr70;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr71;
	} else
		goto tr71;
	goto tr68;
tr69:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st39;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
#line 11338 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr72;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr72;
	goto tr68;
tr72:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
#line 11374 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr52;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr36;
	} else
		goto tr36;
	goto tr73;
tr74:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
#line 11406 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr75;
		case 46: goto st32;
		case 95: goto st25;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr75;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 65 )
			goto st25;
	} else
		goto st41;
	goto tr73;
tr75:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st42;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
#line 11457 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr77;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr78;
	} else if ( (*p) >= 9 )
		goto st24;
	goto tr55;
tr70:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st43;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
#line 11491 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr72;
		case 46: goto st35;
		case 95: goto st26;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr72;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 65 )
			goto st26;
	} else
		goto st43;
	goto tr68;
tr449:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st230;
tr459:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st230;
st230:
	if ( ++p == pe )
		goto _test_eof230;
case 230:
#line 11630 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto tr375;
		case 66: goto tr362;
		case 71: goto tr361;
		case 78: goto tr361;
		case 84: goto tr361;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st99;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr457;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr361;
	goto tr456;
tr457:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st231;
st231:
	if ( ++p == pe )
		goto _test_eof231;
case 231:
#line 11703 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr459;
		case 32: goto tr458;
		case 46: goto st232;
		case 66: goto st44;
		case 71: goto st190;
		case 78: goto st190;
		case 84: goto st190;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr458;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st231;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st190;
	goto tr456;
st232:
	if ( ++p == pe )
		goto _test_eof232;
case 232:
	switch( (*p) ) {
		case 10: goto tr176;
		case 32: goto tr175;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr175;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st233;
	} else
		goto st13;
	goto tr221;
st233:
	if ( ++p == pe )
		goto _test_eof233;
case 233:
	switch( (*p) ) {
		case 10: goto tr463;
		case 32: goto tr462;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr462;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st233;
	} else
		goto st13;
	goto tr221;
tr362:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
#line 11820 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr81;
		case 95: goto st44;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr81;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st44;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st13;
	goto tr80;
tr81:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
#line 11875 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr69;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr84;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr85;
	} else
		goto tr85;
	goto tr83;
tr84:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st46;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
#line 11911 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr86;
		case 46: goto st35;
		case 95: goto st37;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr86;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 65 )
			goto st37;
	} else
		goto st46;
	goto tr83;
tr86:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st47;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
#line 11963 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr69;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr89;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr71;
	} else
		goto tr71;
	goto tr88;
tr89:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 11999 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr90;
		case 46: goto st32;
		case 95: goto st26;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr90;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 65 )
			goto st26;
	} else
		goto st48;
	goto tr88;
tr90:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
#line 12056 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr92;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr93;
		} else if ( (*p) >= 9 )
			goto st24;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr73;
st234:
	if ( ++p == pe )
		goto _test_eof234;
case 234:
	switch( (*p) ) {
		case 10: goto tr176;
		case 32: goto tr175;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr175;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st235;
	} else
		goto st13;
	goto tr224;
st235:
	if ( ++p == pe )
		goto _test_eof235;
case 235:
	switch( (*p) ) {
		case 10: goto tr226;
		case 32: goto tr225;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr225;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st235;
	} else
		goto st13;
	goto tr224;
st236:
	if ( ++p == pe )
		goto _test_eof236;
case 236:
	switch( (*p) ) {
		case 10: goto tr176;
		case 32: goto tr175;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr175;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st237;
	} else
		goto st13;
	goto tr195;
st237:
	if ( ++p == pe )
		goto _test_eof237;
case 237:
	switch( (*p) ) {
		case 10: goto tr198;
		case 32: goto tr197;
		case 46: goto st101;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr197;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st237;
	} else
		goto st13;
	goto tr195;
tr258:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
#line 12208 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr96;
		case 95: goto st50;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr96;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st50;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st13;
	goto tr95;
tr96:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
#line 12257 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr98;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr85;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr85;
	} else
		goto tr85;
	goto tr65;
tr98:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st52;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
#line 12285 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr49;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr49;
	goto tr48;
tr306:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
#line 12314 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr100;
		case 95: goto st53;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr100;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st53;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st13;
	goto tr99;
tr100:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 12363 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr103;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr104;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr105;
	} else
		goto tr105;
	goto tr102;
tr103:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st55;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
#line 12391 "ped_reader.c"
	if ( (*p) == 32 )
		goto tr107;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr107;
	goto tr106;
tr107:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 12421 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr35;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr109;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr36;
	} else
		goto tr36;
	goto tr108;
tr109:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 12449 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr110;
		case 46: goto st59;
		case 95: goto st25;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr110;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 65 )
			goto st25;
	} else
		goto st57;
	goto tr108;
tr110:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
#line 12488 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr77;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr112;
	} else if ( (*p) >= 9 )
		goto st24;
	goto tr15;
tr112:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st238;
st238:
	if ( ++p == pe )
		goto _test_eof238;
case 238:
#line 12514 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto tr188;
		case 46: goto st239;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st238;
	} else if ( (*p) >= 9 )
		goto tr188;
	goto tr15;
st239:
	if ( ++p == pe )
		goto _test_eof239;
case 239:
	switch( (*p) ) {
		case 10: goto tr173;
		case 32: goto st99;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st240;
	} else if ( (*p) >= 9 )
		goto st99;
	goto tr15;
st240:
	if ( ++p == pe )
		goto _test_eof240;
case 240:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto tr188;
		case 46: goto st99;
		case 65: goto st99;
		case 67: goto st99;
		case 71: goto st99;
		case 78: goto st99;
		case 84: goto st99;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st240;
	} else if ( (*p) >= 9 )
		goto tr188;
	goto tr15;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st60;
	goto tr19;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	if ( (*p) == 32 )
		goto tr42;
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st60;
	} else if ( (*p) >= 9 )
		goto tr42;
	goto tr19;
tr104:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
#line 12608 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr114;
		case 46: goto st35;
		case 95: goto st65;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr114;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 65 )
			goto st65;
	} else
		goto st61;
	goto tr102;
tr114:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
#line 12654 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr98;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr118;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr71;
	} else
		goto tr71;
	goto tr117;
tr118:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
#line 12686 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr119;
		case 46: goto st59;
		case 95: goto st26;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr119;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 65 )
			goto st26;
	} else
		goto st63;
	goto tr117;
tr119:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 12731 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr92;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr121;
		} else if ( (*p) >= 9 )
			goto st24;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr51;
tr121:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st241;
st241:
	if ( ++p == pe )
		goto _test_eof241;
case 241:
#line 12770 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr215;
		case 32: goto tr214;
		case 46: goto st239;
		case 66: goto st25;
		case 71: goto st152;
		case 78: goto st152;
		case 84: goto st152;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st241;
		} else if ( (*p) >= 9 )
			goto tr214;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st152;
	goto tr51;
tr105:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 12810 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr123;
		case 95: goto st65;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr123;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 65 )
			goto st65;
	} else
		goto st65;
	goto tr122;
tr123:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
#line 12843 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr98;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr71;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr71;
	} else
		goto tr71;
	goto tr48;
tr334:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st242;
st242:
	if ( ++p == pe )
		goto _test_eof242;
case 242:
#line 12880 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr337;
		case 32: goto tr336;
		case 46: goto st101;
		case 66: goto st67;
		case 71: goto st242;
		case 78: goto st242;
		case 84: goto st242;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr336;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st242;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st242;
	goto tr124;
tr335:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 12941 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr125;
		case 95: goto st67;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr125;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st67;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st13;
	goto tr124;
tr125:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
#line 12992 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr52;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr128;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr129;
	} else
		goto tr129;
	goto tr127;
tr128:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
#line 13024 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr130;
		case 46: goto st35;
		case 95: goto st74;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr130;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 65 )
			goto st74;
	} else
		goto st69;
	goto tr127;
tr130:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
#line 13070 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr103;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr134;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr33;
	} else
		goto tr33;
	goto tr133;
tr134:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
#line 13102 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr135;
		case 46: goto st32;
		case 95: goto st73;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr135;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 65 )
			goto st73;
	} else
		goto st71;
	goto tr133;
tr135:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st72;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
#line 13153 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr138;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr139;
		} else if ( (*p) >= 9 )
			goto st24;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr108;
tr33:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
#line 13188 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr34;
		case 95: goto st73;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr34;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 65 )
			goto st73;
	} else
		goto st73;
	goto tr7;
tr129:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
#line 13219 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr141;
		case 95: goto st74;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr141;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 65 )
			goto st74;
	} else
		goto st74;
	goto tr140;
tr141:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st75;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
#line 13252 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr103;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr142;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr33;
	} else
		goto tr33;
	goto tr106;
tr142:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st76;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
#line 13280 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr107;
		case 46: goto st35;
		case 95: goto st73;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr107;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 65 )
			goto st73;
	} else
		goto st76;
	goto tr106;
tr233:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st243;
st243:
	if ( ++p == pe )
		goto _test_eof243;
case 243:
#line 13317 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr236;
		case 32: goto tr235;
		case 46: goto st101;
		case 66: goto st77;
		case 71: goto st243;
		case 78: goto st243;
		case 84: goto st243;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr235;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st243;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st243;
	goto tr144;
tr234:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st77;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
#line 13374 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr145;
		case 95: goto st77;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr145;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st77;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st13;
	goto tr144;
tr145:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st78;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
#line 13419 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr35;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr129;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr129;
	} else
		goto tr129;
	goto tr140;
tr387:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st244;
st244:
	if ( ++p == pe )
		goto _test_eof244;
case 244:
#line 13460 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr470;
		case 32: goto tr469;
		case 46: goto st232;
		case 66: goto st77;
		case 71: goto st243;
		case 78: goto st243;
		case 84: goto st243;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr469;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st244;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st243;
	goto tr385;
tr376:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st245;
st245:
	if ( ++p == pe )
		goto _test_eof245;
case 245:
#line 13529 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr473;
		case 32: goto tr472;
		case 46: goto st232;
		case 66: goto st67;
		case 71: goto st242;
		case 78: goto st242;
		case 84: goto st242;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr472;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st245;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st242;
	goto tr374;
tr202:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
#line 13586 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr148;
		case 95: goto st79;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr148;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st79;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st13;
	goto tr147;
tr148:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st80;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
#line 13629 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr32;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr105;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr105;
	} else
		goto tr105;
	goto tr122;
tr276:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st246;
st246:
	if ( ++p == pe )
		goto _test_eof246;
case 246:
#line 13662 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr293;
		case 32: goto tr292;
		case 46: goto st101;
		case 66: goto st81;
		case 71: goto st246;
		case 78: goto st246;
		case 84: goto st246;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr292;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st246;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st246;
	goto tr150;
tr277:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st81;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
#line 13719 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr151;
		case 95: goto st81;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr151;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st81;
		} else
			goto st13;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st13;
	goto tr150;
tr151:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
#line 13764 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr38;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr154;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr153;
tr154:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
#line 13792 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr155;
		case 46: goto st35;
		case 95: goto st15;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr155;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 65 )
			goto st15;
	} else
		goto st83;
	goto tr153;
tr155:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st84;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
#line 13832 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr32;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr158;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr33;
	} else
		goto tr33;
	goto tr157;
tr158:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
#line 13860 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr159;
		case 46: goto st59;
		case 95: goto st73;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr159;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 65 )
			goto st73;
	} else
		goto st85;
	goto tr157;
tr159:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st86;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
#line 13899 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto st24;
		case 46: goto tr138;
		case 66: goto tr36;
		case 71: goto tr94;
		case 78: goto tr94;
		case 84: goto tr94;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr94;
		} else if ( (*p) >= 9 )
			goto st24;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr94;
	goto tr11;
tr183:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st247;
st247:
	if ( ++p == pe )
		goto _test_eof247;
case 247:
#line 13934 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr185;
		case 32: goto tr184;
		case 46: goto st99;
		case 66: goto st73;
		case 71: goto st247;
		case 78: goto st247;
		case 84: goto st247;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st247;
		} else if ( (*p) >= 9 )
			goto tr184;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st247;
	goto tr7;
tr275:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st248;
st248:
	if ( ++p == pe )
		goto _test_eof248;
case 248:
#line 13987 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr476;
		case 32: goto tr475;
		case 46: goto st232;
		case 66: goto st81;
		case 71: goto st246;
		case 78: goto st246;
		case 84: goto st246;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr475;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st248;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st246;
	goto tr273;
tr267:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st249;
st249:
	if ( ++p == pe )
		goto _test_eof249;
case 249:
#line 14043 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr269;
		case 32: goto tr268;
		case 46: goto st239;
		case 66: goto st26;
		case 71: goto st195;
		case 78: goto st195;
		case 84: goto st195;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st249;
		} else if ( (*p) >= 9 )
			goto tr268;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st195;
	goto tr68;
tr242:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st250;
st250:
	if ( ++p == pe )
		goto _test_eof250;
case 250:
#line 14083 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr244;
		case 32: goto tr243;
		case 46: goto st239;
		case 66: goto st73;
		case 71: goto st247;
		case 78: goto st247;
		case 84: goto st247;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st250;
		} else if ( (*p) >= 9 )
			goto tr243;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st247;
	goto tr106;
tr223:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st251;
st251:
	if ( ++p == pe )
		goto _test_eof251;
case 251:
#line 14132 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr463;
		case 32: goto tr462;
		case 46: goto st232;
		case 65: goto st101;
		case 67: goto st101;
		case 71: goto st101;
		case 78: goto st101;
		case 84: goto st101;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr462;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st251;
	} else
		goto st13;
	goto tr221;
tr39:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
#line 14165 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr40;
		case 46: goto st35;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st87;
	} else if ( (*p) >= 9 )
		goto tr40;
	goto tr15;
tr172:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
	goto st252;
st252:
	if ( ++p == pe )
		goto _test_eof252;
case 252:
#line 14191 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr482;
		case 32: goto tr481;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr481;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr0;
tr481:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st253;
tr507:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st253;
st253:
	if ( ++p == pe )
		goto _test_eof253;
case 253:
#line 14239 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 66: goto tr29;
		case 71: goto tr484;
		case 78: goto tr484;
		case 84: goto tr484;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr484;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr484;
	goto tr3;
tr484:
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st254;
st254:
	if ( ++p == pe )
		goto _test_eof254;
case 254:
#line 14275 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr486;
		case 32: goto tr485;
		case 46: goto st97;
		case 66: goto st15;
		case 71: goto st254;
		case 78: goto st254;
		case 84: goto st254;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st254;
		} else if ( (*p) >= 9 )
			goto tr485;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st254;
	goto tr3;
tr485:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st255;
tr538:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st255;
st255:
	if ( ++p == pe )
		goto _test_eof255;
case 255:
#line 14326 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr488;
		case 66: goto tr33;
		case 71: goto tr489;
		case 78: goto tr489;
		case 84: goto tr489;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr489;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr489;
	goto tr7;
tr488:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st256;
st256:
	if ( ++p == pe )
		goto _test_eof256;
case 256:
#line 14362 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr491;
		case 32: goto tr490;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr490;
	goto tr7;
tr490:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st257;
tr591:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st257;
st257:
	if ( ++p == pe )
		goto _test_eof257;
case 257:
#line 14408 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr492;
		case 66: goto tr36;
		case 71: goto tr493;
		case 78: goto tr493;
		case 84: goto tr493;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr493;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr493;
	goto tr11;
tr492:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st258;
st258:
	if ( ++p == pe )
		goto _test_eof258;
case 258:
#line 14444 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr495;
		case 32: goto tr494;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr494;
	goto tr11;
tr494:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st259;
tr554:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st259;
st259:
	if ( ++p == pe )
		goto _test_eof259;
case 259:
#line 14490 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr496;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr497;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr15;
tr496:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st260;
st260:
	if ( ++p == pe )
		goto _test_eof260;
case 260:
#line 14517 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr499;
		case 32: goto tr498;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr498;
	goto tr15;
tr498:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st261;
tr528:
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st261;
st261:
	if ( ++p == pe )
		goto _test_eof261;
case 261:
#line 14575 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr500;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr19;
tr500:
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st262;
st262:
	if ( ++p == pe )
		goto _test_eof262;
case 262:
#line 14602 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr502;
		case 32: goto tr501;
		case 46: goto st263;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st262;
	} else if ( (*p) >= 9 )
		goto tr501;
	goto tr19;
st263:
	if ( ++p == pe )
		goto _test_eof263;
case 263:
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st264;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr19;
st264:
	if ( ++p == pe )
		goto _test_eof264;
case 264:
	switch( (*p) ) {
		case 10: goto tr502;
		case 32: goto tr501;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st264;
	} else if ( (*p) >= 9 )
		goto tr501;
	goto tr19;
tr499:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st265;
tr529:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st265;
st265:
	if ( ++p == pe )
		goto _test_eof265;
case 265:
#line 14740 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr172;
		case 65: goto tr172;
		case 67: goto tr172;
		case 71: goto tr172;
		case 78: goto tr172;
		case 84: goto tr172;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st97;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr506;
	} else
		goto tr171;
	goto tr195;
tr506:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st266;
st266:
	if ( ++p == pe )
		goto _test_eof266;
case 266:
#line 14782 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr508;
		case 32: goto tr507;
		case 46: goto st387;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr507;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st266;
	} else
		goto st13;
	goto tr195;
tr482:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st267;
tr508:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st267;
st267:
	if ( ++p == pe )
		goto _test_eof267;
case 267:
#line 14870 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr172;
		case 66: goto tr202;
		case 71: goto tr511;
		case 78: goto tr511;
		case 84: goto tr511;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr511;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr511;
	goto tr147;
tr511:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
	goto st268;
st268:
	if ( ++p == pe )
		goto _test_eof268;
case 268:
#line 14927 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr513;
		case 32: goto tr512;
		case 46: goto st252;
		case 66: goto st79;
		case 71: goto st268;
		case 78: goto st268;
		case 84: goto st268;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr512;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st268;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st268;
	goto tr147;
tr512:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st269;
tr562:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st269;
st269:
	if ( ++p == pe )
		goto _test_eof269;
case 269:
#line 14998 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr488;
		case 66: goto tr105;
		case 71: goto tr515;
		case 78: goto tr515;
		case 84: goto tr515;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr515;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr515;
	goto tr122;
tr515:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st270;
st270:
	if ( ++p == pe )
		goto _test_eof270;
case 270:
#line 15038 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr517;
		case 32: goto tr516;
		case 46: goto st97;
		case 66: goto st65;
		case 71: goto st270;
		case 78: goto st270;
		case 84: goto st270;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st270;
		} else if ( (*p) >= 9 )
			goto tr516;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st270;
	goto tr122;
tr516:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st271;
tr630:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st271;
st271:
	if ( ++p == pe )
		goto _test_eof271;
case 271:
#line 15101 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr519;
		case 66: goto tr71;
		case 71: goto tr520;
		case 78: goto tr520;
		case 84: goto tr520;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr520;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr520;
	goto tr48;
tr519:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st272;
st272:
	if ( ++p == pe )
		goto _test_eof272;
case 272:
#line 15141 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr522;
		case 32: goto tr521;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr521;
	goto tr48;
tr521:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st273;
tr617:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st273;
st273:
	if ( ++p == pe )
		goto _test_eof273;
case 273:
#line 15199 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr523;
		case 66: goto tr36;
		case 71: goto tr493;
		case 78: goto tr493;
		case 84: goto tr493;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr524;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr493;
	goto tr51;
tr523:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st274;
st274:
	if ( ++p == pe )
		goto _test_eof274;
case 274:
#line 15239 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr526;
		case 32: goto tr525;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr525;
	goto tr51;
tr525:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st275;
tr578:
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st275;
st275:
	if ( ++p == pe )
		goto _test_eof275;
case 275:
#line 15309 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr496;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr527;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr55;
tr527:
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st276;
st276:
	if ( ++p == pe )
		goto _test_eof276;
case 276:
#line 15340 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr529;
		case 32: goto tr528;
		case 46: goto st277;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st276;
	} else if ( (*p) >= 9 )
		goto tr528;
	goto tr55;
st277:
	if ( ++p == pe )
		goto _test_eof277;
case 277:
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st278;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr55;
st278:
	if ( ++p == pe )
		goto _test_eof278;
case 278:
	switch( (*p) ) {
		case 10: goto tr529;
		case 32: goto tr528;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st278;
	} else if ( (*p) >= 9 )
		goto tr528;
	goto tr55;
tr526:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st279;
tr579:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st279;
st279:
	if ( ++p == pe )
		goto _test_eof279;
case 279:
#line 15490 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr533;
		case 65: goto tr172;
		case 67: goto tr172;
		case 71: goto tr172;
		case 78: goto tr172;
		case 84: goto tr172;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st97;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr534;
	} else
		goto tr171;
	goto tr221;
tr533:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st280;
st280:
	if ( ++p == pe )
		goto _test_eof280;
case 280:
#line 15532 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr536;
		case 32: goto tr535;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr535;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr224;
tr535:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st281;
tr739:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st281;
st281:
	if ( ++p == pe )
		goto _test_eof281;
case 281:
#line 15604 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 66: goto tr29;
		case 71: goto tr484;
		case 78: goto tr484;
		case 84: goto tr484;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr537;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr484;
	goto tr227;
tr537:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st282;
st282:
	if ( ++p == pe )
		goto _test_eof282;
case 282:
#line 15644 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr539;
		case 32: goto tr538;
		case 46: goto st263;
		case 66: goto st15;
		case 71: goto st254;
		case 78: goto st254;
		case 84: goto st254;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st282;
		} else if ( (*p) >= 9 )
			goto tr538;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st254;
	goto tr227;
tr486:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st283;
tr539:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st283;
st283:
	if ( ++p == pe )
		goto _test_eof283;
case 283:
#line 15735 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr541;
		case 66: goto tr234;
		case 71: goto tr542;
		case 78: goto tr542;
		case 84: goto tr542;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr542;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr542;
	goto tr144;
tr541:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st284;
st284:
	if ( ++p == pe )
		goto _test_eof284;
case 284:
#line 15792 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr544;
		case 32: goto tr543;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr543;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr144;
tr543:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st285;
tr604:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st285;
st285:
	if ( ++p == pe )
		goto _test_eof285;
case 285:
#line 15852 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr492;
		case 66: goto tr129;
		case 71: goto tr545;
		case 78: goto tr545;
		case 84: goto tr545;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr545;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr545;
	goto tr140;
tr545:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st286;
st286:
	if ( ++p == pe )
		goto _test_eof286;
case 286:
#line 15892 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr547;
		case 32: goto tr546;
		case 46: goto st97;
		case 66: goto st74;
		case 71: goto st286;
		case 78: goto st286;
		case 84: goto st286;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st286;
		} else if ( (*p) >= 9 )
			goto tr546;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st286;
	goto tr140;
tr546:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st287;
tr682:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st287;
st287:
	if ( ++p == pe )
		goto _test_eof287;
case 287:
#line 15955 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr549;
		case 66: goto tr33;
		case 71: goto tr489;
		case 78: goto tr489;
		case 84: goto tr489;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr550;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr489;
	goto tr106;
tr549:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st288;
st288:
	if ( ++p == pe )
		goto _test_eof288;
case 288:
#line 15995 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr552;
		case 32: goto tr551;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr551;
	goto tr106;
tr551:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st289;
tr642:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st289;
st289:
	if ( ++p == pe )
		goto _test_eof289;
case 289:
#line 16065 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr492;
		case 66: goto tr36;
		case 71: goto tr493;
		case 78: goto tr493;
		case 84: goto tr493;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr553;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr493;
	goto tr108;
tr553:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st290;
st290:
	if ( ++p == pe )
		goto _test_eof290;
case 290:
#line 16105 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr555;
		case 32: goto tr554;
		case 46: goto st263;
		case 66: goto st25;
		case 71: goto st303;
		case 78: goto st303;
		case 84: goto st303;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st290;
		} else if ( (*p) >= 9 )
			goto tr554;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st303;
	goto tr108;
tr495:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st291;
tr555:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st291;
st291:
	if ( ++p == pe )
		goto _test_eof291;
case 291:
#line 16200 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr533;
		case 65: goto tr172;
		case 67: goto tr172;
		case 71: goto tr172;
		case 78: goto tr172;
		case 84: goto tr172;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st97;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto tr171;
		} else if ( (*p) >= 48 )
			goto tr558;
	} else
		goto tr171;
	goto tr224;
tr558:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st292;
st292:
	if ( ++p == pe )
		goto _test_eof292;
case 292:
#line 16242 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr536;
		case 32: goto tr535;
		case 46: goto st385;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr535;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st292;
	} else
		goto st13;
	goto tr224;
tr536:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st293;
tr740:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st293;
st293:
	if ( ++p == pe )
		goto _test_eof293;
case 293:
#line 16354 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr172;
		case 66: goto tr202;
		case 71: goto tr511;
		case 78: goto tr511;
		case 84: goto tr511;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr561;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr511;
	goto tr252;
tr561:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st294;
st294:
	if ( ++p == pe )
		goto _test_eof294;
case 294:
#line 16415 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr563;
		case 32: goto tr562;
		case 46: goto st387;
		case 66: goto st79;
		case 71: goto st268;
		case 78: goto st268;
		case 84: goto st268;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr562;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st294;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st268;
	goto tr252;
tr513:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st295;
tr563:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st295;
st295:
	if ( ++p == pe )
		goto _test_eof295;
case 295:
#line 16526 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr541;
		case 66: goto tr258;
		case 71: goto tr565;
		case 78: goto tr565;
		case 84: goto tr565;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr565;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr565;
	goto tr95;
tr565:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st296;
st296:
	if ( ++p == pe )
		goto _test_eof296;
case 296:
#line 16587 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr567;
		case 32: goto tr566;
		case 46: goto st252;
		case 66: goto st50;
		case 71: goto st296;
		case 78: goto st296;
		case 84: goto st296;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr566;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st296;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st296;
	goto tr95;
tr566:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st297;
tr654:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st297;
st297:
	if ( ++p == pe )
		goto _test_eof297;
case 297:
#line 16670 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr519;
		case 66: goto tr85;
		case 71: goto tr569;
		case 78: goto tr569;
		case 84: goto tr569;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr569;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr569;
	goto tr65;
tr569:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st298;
st298:
	if ( ++p == pe )
		goto _test_eof298;
case 298:
#line 16714 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr571;
		case 32: goto tr570;
		case 46: goto st97;
		case 66: goto st37;
		case 71: goto st298;
		case 78: goto st298;
		case 84: goto st298;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st298;
		} else if ( (*p) >= 9 )
			goto tr570;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st298;
	goto tr65;
tr570:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st299;
tr696:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st299;
st299:
	if ( ++p == pe )
		goto _test_eof299;
case 299:
#line 16789 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr573;
		case 66: goto tr71;
		case 71: goto tr520;
		case 78: goto tr520;
		case 84: goto tr520;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr574;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr520;
	goto tr68;
tr573:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st300;
st300:
	if ( ++p == pe )
		goto _test_eof300;
case 300:
#line 16833 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr576;
		case 32: goto tr575;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 9 )
		goto tr575;
	goto tr68;
tr575:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st301;
tr666:
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st301;
st301:
	if ( ++p == pe )
		goto _test_eof301;
case 301:
#line 16915 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr523;
		case 66: goto tr36;
		case 71: goto tr493;
		case 78: goto tr493;
		case 84: goto tr493;
		case 95: goto tr36;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr577;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr36;
		} else if ( (*p) >= 68 )
			goto tr36;
	} else
		goto tr493;
	goto tr73;
tr577:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st302;
st302:
	if ( ++p == pe )
		goto _test_eof302;
case 302:
#line 16959 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr579;
		case 32: goto tr578;
		case 46: goto st277;
		case 66: goto st25;
		case 71: goto st303;
		case 78: goto st303;
		case 84: goto st303;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st302;
		} else if ( (*p) >= 9 )
			goto tr578;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st303;
	goto tr73;
tr493:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st303;
st303:
	if ( ++p == pe )
		goto _test_eof303;
case 303:
#line 16995 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr495;
		case 32: goto tr494;
		case 46: goto st97;
		case 66: goto st25;
		case 71: goto st303;
		case 78: goto st303;
		case 84: goto st303;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st303;
		} else if ( (*p) >= 9 )
			goto tr494;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st303;
	goto tr11;
tr576:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st304;
tr667:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st304;
st304:
	if ( ++p == pe )
		goto _test_eof304;
case 304:
#line 17126 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr581;
		case 66: goto tr277;
		case 71: goto tr583;
		case 78: goto tr583;
		case 84: goto tr583;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr582;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr583;
	goto tr273;
tr581:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st305;
st305:
	if ( ++p == pe )
		goto _test_eof305;
case 305:
#line 17187 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr585;
		case 32: goto tr584;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr584;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr278;
tr584:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st306;
tr750:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st306;
st306:
	if ( ++p == pe )
		goto _test_eof306;
case 306:
#line 17271 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr496;
		case 66: goto tr29;
		case 71: goto tr484;
		case 78: goto tr484;
		case 84: goto tr484;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr586;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr484;
	goto tr281;
tr586:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st307;
st307:
	if ( ++p == pe )
		goto _test_eof307;
case 307:
#line 17315 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr588;
		case 32: goto tr587;
		case 46: goto st277;
		case 66: goto st15;
		case 71: goto st254;
		case 78: goto st254;
		case 84: goto st254;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st307;
		} else if ( (*p) >= 9 )
			goto tr587;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st254;
	goto tr281;
tr599:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st308;
tr587:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st308;
st308:
	if ( ++p == pe )
		goto _test_eof308;
case 308:
#line 17390 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr488;
		case 66: goto tr33;
		case 71: goto tr489;
		case 78: goto tr489;
		case 84: goto tr489;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr590;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr489;
	goto tr157;
tr590:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st309;
st309:
	if ( ++p == pe )
		goto _test_eof309;
case 309:
#line 17430 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr592;
		case 32: goto tr591;
		case 46: goto st263;
		case 66: goto st73;
		case 71: goto st396;
		case 78: goto st396;
		case 84: goto st396;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st309;
		} else if ( (*p) >= 9 )
			goto tr591;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st396;
	goto tr157;
tr491:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st310;
tr592:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st310;
st310:
	if ( ++p == pe )
		goto _test_eof310;
case 310:
#line 17525 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr595;
		case 66: goto tr277;
		case 71: goto tr583;
		case 78: goto tr583;
		case 84: goto tr583;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr583;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr583;
	goto tr150;
tr595:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st311;
st311:
	if ( ++p == pe )
		goto _test_eof311;
case 311:
#line 17582 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr597;
		case 32: goto tr596;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr596;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr150;
tr596:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st312;
tr646:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st312;
st312:
	if ( ++p == pe )
		goto _test_eof312;
case 312:
#line 17642 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr496;
		case 66: goto tr29;
		case 71: goto tr484;
		case 78: goto tr484;
		case 84: goto tr484;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr598;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr29;
		} else if ( (*p) >= 68 )
			goto tr29;
	} else
		goto tr484;
	goto tr153;
tr598:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st313;
st313:
	if ( ++p == pe )
		goto _test_eof313;
case 313:
#line 17682 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr600;
		case 32: goto tr599;
		case 46: goto st393;
		case 66: goto st15;
		case 71: goto st254;
		case 78: goto st254;
		case 84: goto st254;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st313;
		} else if ( (*p) >= 9 )
			goto tr599;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 68 )
			goto st15;
	} else
		goto st254;
	goto tr153;
tr588:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st314;
tr600:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st314;
st314:
	if ( ++p == pe )
		goto _test_eof314;
case 314:
#line 17797 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr541;
		case 66: goto tr234;
		case 71: goto tr542;
		case 78: goto tr542;
		case 84: goto tr542;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr603;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr542;
	goto tr299;
tr603:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st315;
st315:
	if ( ++p == pe )
		goto _test_eof315;
case 315:
#line 17858 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr605;
		case 32: goto tr604;
		case 46: goto st387;
		case 66: goto st77;
		case 71: goto st390;
		case 78: goto st390;
		case 84: goto st390;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr604;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st315;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st390;
	goto tr299;
tr544:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st316;
tr605:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st316;
st316:
	if ( ++p == pe )
		goto _test_eof316;
case 316:
#line 17973 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr595;
		case 66: goto tr306;
		case 71: goto tr608;
		case 78: goto tr608;
		case 84: goto tr608;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr608;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr608;
	goto tr99;
tr608:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st317;
st317:
	if ( ++p == pe )
		goto _test_eof317;
case 317:
#line 18034 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr610;
		case 32: goto tr609;
		case 46: goto st252;
		case 66: goto st53;
		case 71: goto st317;
		case 78: goto st317;
		case 84: goto st317;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr609;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st317;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st317;
	goto tr99;
tr609:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st318;
tr688:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st318;
st318:
	if ( ++p == pe )
		goto _test_eof318;
case 318:
#line 18117 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr549;
		case 66: goto tr105;
		case 71: goto tr515;
		case 78: goto tr515;
		case 84: goto tr515;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr612;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr515;
	goto tr102;
tr612:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st319;
st319:
	if ( ++p == pe )
		goto _test_eof319;
case 319:
#line 18161 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr614;
		case 32: goto tr613;
		case 46: goto st393;
		case 66: goto st65;
		case 71: goto st270;
		case 78: goto st270;
		case 84: goto st270;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st319;
		} else if ( (*p) >= 9 )
			goto tr613;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st270;
	goto tr102;
tr613:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st320;
tr707:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st320;
st320:
	if ( ++p == pe )
		goto _test_eof320;
case 320:
#line 18248 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr519;
		case 66: goto tr71;
		case 71: goto tr520;
		case 78: goto tr520;
		case 84: goto tr520;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr616;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr520;
	goto tr117;
tr616:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st321;
st321:
	if ( ++p == pe )
		goto _test_eof321;
case 321:
#line 18292 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr618;
		case 32: goto tr617;
		case 46: goto st263;
		case 66: goto st26;
		case 71: goto st346;
		case 78: goto st346;
		case 84: goto st346;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st321;
		} else if ( (*p) >= 9 )
			goto tr617;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st346;
	goto tr117;
tr522:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st322;
tr618:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st322;
st322:
	if ( ++p == pe )
		goto _test_eof322;
case 322:
#line 18399 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr581;
		case 66: goto tr277;
		case 71: goto tr583;
		case 78: goto tr583;
		case 84: goto tr583;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr621;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr583;
	goto tr278;
tr621:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st323;
st323:
	if ( ++p == pe )
		goto _test_eof323;
case 323:
#line 18460 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr585;
		case 32: goto tr584;
		case 46: goto st385;
		case 66: goto st81;
		case 71: goto st395;
		case 78: goto st395;
		case 84: goto st395;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr584;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st323;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st395;
	goto tr278;
tr585:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st324;
tr751:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st324;
st324:
	if ( ++p == pe )
		goto _test_eof324;
case 324:
#line 18599 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr533;
		case 66: goto tr202;
		case 71: goto tr511;
		case 78: goto tr511;
		case 84: goto tr511;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr624;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr511;
	goto tr322;
tr624:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st325;
st325:
	if ( ++p == pe )
		goto _test_eof325;
case 325:
#line 18664 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr626;
		case 32: goto tr625;
		case 46: goto st383;
		case 66: goto st79;
		case 71: goto st268;
		case 78: goto st268;
		case 84: goto st268;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr625;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st325;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st268;
	goto tr322;
tr625:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st326;
tr650:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st326;
st326:
	if ( ++p == pe )
		goto _test_eof326;
case 326:
#line 18759 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr488;
		case 66: goto tr105;
		case 71: goto tr515;
		case 78: goto tr515;
		case 84: goto tr515;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr629;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr515;
	goto tr328;
tr629:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st327;
st327:
	if ( ++p == pe )
		goto _test_eof327;
case 327:
#line 18803 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr631;
		case 32: goto tr630;
		case 46: goto st263;
		case 66: goto st65;
		case 71: goto st270;
		case 78: goto st270;
		case 84: goto st270;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st327;
		} else if ( (*p) >= 9 )
			goto tr630;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st270;
	goto tr328;
tr517:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st328;
tr631:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st328;
st328:
	if ( ++p == pe )
		goto _test_eof328;
case 328:
#line 18906 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr633;
		case 66: goto tr335;
		case 71: goto tr634;
		case 78: goto tr634;
		case 84: goto tr634;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr634;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr634;
	goto tr124;
tr633:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st329;
st329:
	if ( ++p == pe )
		goto _test_eof329;
case 329:
#line 18967 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr636;
		case 32: goto tr635;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr635;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr124;
tr635:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st330;
tr711:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st330;
st330:
	if ( ++p == pe )
		goto _test_eof330;
case 330:
#line 19039 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr523;
		case 66: goto tr129;
		case 71: goto tr545;
		case 78: goto tr545;
		case 84: goto tr545;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr637;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr545;
	goto tr127;
tr637:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st331;
st331:
	if ( ++p == pe )
		goto _test_eof331;
case 331:
#line 19083 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr639;
		case 32: goto tr638;
		case 46: goto st393;
		case 66: goto st74;
		case 71: goto st286;
		case 78: goto st286;
		case 84: goto st286;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st331;
		} else if ( (*p) >= 9 )
			goto tr638;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st286;
	goto tr127;
tr638:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st332;
tr674:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st332;
st332:
	if ( ++p == pe )
		goto _test_eof332;
case 332:
#line 19170 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr549;
		case 66: goto tr33;
		case 71: goto tr489;
		case 78: goto tr489;
		case 84: goto tr489;
		case 95: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr641;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr33;
		} else if ( (*p) >= 68 )
			goto tr33;
	} else
		goto tr489;
	goto tr133;
tr641:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st333;
st333:
	if ( ++p == pe )
		goto _test_eof333;
case 333:
#line 19214 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr643;
		case 32: goto tr642;
		case 46: goto st277;
		case 66: goto st73;
		case 71: goto st396;
		case 78: goto st396;
		case 84: goto st396;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st333;
		} else if ( (*p) >= 9 )
			goto tr642;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st396;
	goto tr133;
tr552:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st334;
tr643:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st334;
st334:
	if ( ++p == pe )
		goto _test_eof334;
case 334:
#line 19333 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr595;
		case 66: goto tr277;
		case 71: goto tr583;
		case 78: goto tr583;
		case 84: goto tr583;
		case 95: goto tr277;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr645;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr277;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr277;
		} else
			goto tr171;
	} else
		goto tr583;
	goto tr346;
tr645:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st335;
st335:
	if ( ++p == pe )
		goto _test_eof335;
case 335:
#line 19394 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr647;
		case 32: goto tr646;
		case 46: goto st387;
		case 66: goto st81;
		case 71: goto st395;
		case 78: goto st395;
		case 84: goto st395;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr646;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st335;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st395;
	goto tr346;
tr597:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st336;
tr647:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st336;
st336:
	if ( ++p == pe )
		goto _test_eof336;
case 336:
#line 19509 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr533;
		case 66: goto tr202;
		case 71: goto tr511;
		case 78: goto tr511;
		case 84: goto tr511;
		case 95: goto tr202;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr649;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr202;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr202;
		} else
			goto tr171;
	} else
		goto tr511;
	goto tr351;
tr649:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st337;
st337:
	if ( ++p == pe )
		goto _test_eof337;
case 337:
#line 19570 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr651;
		case 32: goto tr650;
		case 46: goto st385;
		case 66: goto st79;
		case 71: goto st268;
		case 78: goto st268;
		case 84: goto st268;
		case 95: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr650;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st337;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st79;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st79;
		} else
			goto st13;
	} else
		goto st268;
	goto tr351;
tr626:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st338;
tr651:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st338;
st338:
	if ( ++p == pe )
		goto _test_eof338;
case 338:
#line 19705 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr541;
		case 66: goto tr258;
		case 71: goto tr565;
		case 78: goto tr565;
		case 84: goto tr565;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr653;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr565;
	goto tr356;
tr653:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st339;
st339:
	if ( ++p == pe )
		goto _test_eof339;
case 339:
#line 19770 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr655;
		case 32: goto tr654;
		case 46: goto st387;
		case 66: goto st50;
		case 71: goto st296;
		case 78: goto st296;
		case 84: goto st296;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr654;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st339;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st296;
	goto tr356;
tr567:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st340;
tr655:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st340;
st340:
	if ( ++p == pe )
		goto _test_eof340;
case 340:
#line 19893 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr633;
		case 66: goto tr362;
		case 71: goto tr657;
		case 78: goto tr657;
		case 84: goto tr657;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr657;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr657;
	goto tr80;
tr657:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st341;
st341:
	if ( ++p == pe )
		goto _test_eof341;
case 341:
#line 19958 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr659;
		case 32: goto tr658;
		case 46: goto st252;
		case 66: goto st44;
		case 71: goto st341;
		case 78: goto st341;
		case 84: goto st341;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr658;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st341;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st341;
	goto tr80;
tr658:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st342;
tr723:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st342;
st342:
	if ( ++p == pe )
		goto _test_eof342;
case 342:
#line 20053 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr573;
		case 66: goto tr85;
		case 71: goto tr569;
		case 78: goto tr569;
		case 84: goto tr569;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr661;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr569;
	goto tr83;
tr661:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st343;
st343:
	if ( ++p == pe )
		goto _test_eof343;
case 343:
#line 20101 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr663;
		case 32: goto tr662;
		case 46: goto st393;
		case 66: goto st37;
		case 71: goto st298;
		case 78: goto st298;
		case 84: goto st298;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st343;
		} else if ( (*p) >= 9 )
			goto tr662;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st298;
	goto tr83;
tr662:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st344;
tr731:
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st344;
st344:
	if ( ++p == pe )
		goto _test_eof344;
case 344:
#line 20200 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr573;
		case 66: goto tr71;
		case 71: goto tr520;
		case 78: goto tr520;
		case 84: goto tr520;
		case 95: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr665;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr71;
		} else if ( (*p) >= 68 )
			goto tr71;
	} else
		goto tr520;
	goto tr88;
tr665:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st345;
st345:
	if ( ++p == pe )
		goto _test_eof345;
case 345:
#line 20248 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr667;
		case 32: goto tr666;
		case 46: goto st277;
		case 66: goto st26;
		case 71: goto st346;
		case 78: goto st346;
		case 84: goto st346;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st345;
		} else if ( (*p) >= 9 )
			goto tr666;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st346;
	goto tr88;
tr520:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st346;
st346:
	if ( ++p == pe )
		goto _test_eof346;
case 346:
#line 20288 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr522;
		case 32: goto tr521;
		case 46: goto st97;
		case 66: goto st26;
		case 71: goto st346;
		case 78: goto st346;
		case 84: goto st346;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st346;
		} else if ( (*p) >= 9 )
			goto tr521;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st346;
	goto tr48;
tr663:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st347;
tr732:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st347;
st347:
	if ( ++p == pe )
		goto _test_eof347;
case 347:
#line 20427 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr669;
		case 66: goto tr335;
		case 71: goto tr634;
		case 78: goto tr634;
		case 84: goto tr634;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr670;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr634;
	goto tr374;
tr669:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st348;
st348:
	if ( ++p == pe )
		goto _test_eof348;
case 348:
#line 20492 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr672;
		case 32: goto tr671;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr671;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr377;
tr671:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st349;
tr746:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st349;
st349:
	if ( ++p == pe )
		goto _test_eof349;
case 349:
#line 20588 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr523;
		case 66: goto tr129;
		case 71: goto tr545;
		case 78: goto tr545;
		case 84: goto tr545;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr673;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr545;
	goto tr380;
tr673:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st350;
st350:
	if ( ++p == pe )
		goto _test_eof350;
case 350:
#line 20636 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr675;
		case 32: goto tr674;
		case 46: goto st277;
		case 66: goto st74;
		case 71: goto st286;
		case 78: goto st286;
		case 84: goto st286;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st350;
		} else if ( (*p) >= 9 )
			goto tr674;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st286;
	goto tr380;
tr639:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st351;
tr675:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st351;
st351:
	if ( ++p == pe )
		goto _test_eof351;
case 351:
#line 20763 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr677;
		case 66: goto tr234;
		case 71: goto tr542;
		case 78: goto tr542;
		case 84: goto tr542;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr678;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr542;
	goto tr385;
tr677:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st352;
st352:
	if ( ++p == pe )
		goto _test_eof352;
case 352:
#line 20824 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr680;
		case 32: goto tr679;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr679;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st252;
	} else
		goto st13;
	goto tr388;
tr679:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st353;
tr743:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st353;
st353:
	if ( ++p == pe )
		goto _test_eof353;
case 353:
#line 20908 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr492;
		case 66: goto tr129;
		case 71: goto tr545;
		case 78: goto tr545;
		case 84: goto tr545;
		case 95: goto tr129;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr681;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr129;
		} else if ( (*p) >= 68 )
			goto tr129;
	} else
		goto tr545;
	goto tr391;
tr681:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st354;
st354:
	if ( ++p == pe )
		goto _test_eof354;
case 354:
#line 20952 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr683;
		case 32: goto tr682;
		case 46: goto st263;
		case 66: goto st74;
		case 71: goto st286;
		case 78: goto st286;
		case 84: goto st286;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st354;
		} else if ( (*p) >= 9 )
			goto tr682;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st74;
		} else if ( (*p) >= 68 )
			goto st74;
	} else
		goto st286;
	goto tr391;
tr547:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st355;
tr683:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st355;
st355:
	if ( ++p == pe )
		goto _test_eof355;
case 355:
#line 21055 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr677;
		case 66: goto tr234;
		case 71: goto tr542;
		case 78: goto tr542;
		case 84: goto tr542;
		case 95: goto tr234;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr685;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr234;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr234;
		} else
			goto tr171;
	} else
		goto tr542;
	goto tr388;
tr685:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st356;
st356:
	if ( ++p == pe )
		goto _test_eof356;
case 356:
#line 21116 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr680;
		case 32: goto tr679;
		case 46: goto st385;
		case 66: goto st77;
		case 71: goto st390;
		case 78: goto st390;
		case 84: goto st390;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr679;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st356;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st390;
	goto tr388;
tr680:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st357;
tr744:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st357;
st357:
	if ( ++p == pe )
		goto _test_eof357;
case 357:
#line 21255 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr595;
		case 66: goto tr306;
		case 71: goto tr608;
		case 78: goto tr608;
		case 84: goto tr608;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr687;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr608;
	goto tr398;
tr687:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st358;
st358:
	if ( ++p == pe )
		goto _test_eof358;
case 358:
#line 21320 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr689;
		case 32: goto tr688;
		case 46: goto st387;
		case 66: goto st53;
		case 71: goto st317;
		case 78: goto st317;
		case 84: goto st317;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr688;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st358;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st317;
	goto tr398;
tr610:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st359;
tr689:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st359;
st359:
	if ( ++p == pe )
		goto _test_eof359;
case 359:
#line 21443 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr677;
		case 66: goto tr258;
		case 71: goto tr565;
		case 78: goto tr565;
		case 84: goto tr565;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr691;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr565;
	goto tr403;
tr691:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st360;
st360:
	if ( ++p == pe )
		goto _test_eof360;
case 360:
#line 21508 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr693;
		case 32: goto tr692;
		case 46: goto st385;
		case 66: goto st50;
		case 71: goto st296;
		case 78: goto st296;
		case 84: goto st296;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr692;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st360;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st296;
	goto tr403;
tr692:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st361;
tr719:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st361;
st361:
	if ( ++p == pe )
		goto _test_eof361;
case 361:
#line 21615 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr519;
		case 66: goto tr85;
		case 71: goto tr569;
		case 78: goto tr569;
		case 84: goto tr569;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr695;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr569;
	goto tr408;
tr695:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st362;
st362:
	if ( ++p == pe )
		goto _test_eof362;
case 362:
#line 21663 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr697;
		case 32: goto tr696;
		case 46: goto st263;
		case 66: goto st37;
		case 71: goto st298;
		case 78: goto st298;
		case 84: goto st298;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st362;
		} else if ( (*p) >= 9 )
			goto tr696;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st298;
	goto tr408;
tr571:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st363;
tr697:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st363;
st363:
	if ( ++p == pe )
		goto _test_eof363;
case 363:
#line 21778 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr669;
		case 66: goto tr335;
		case 71: goto tr634;
		case 78: goto tr634;
		case 84: goto tr634;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr699;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr634;
	goto tr377;
tr699:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st364;
st364:
	if ( ++p == pe )
		goto _test_eof364;
case 364:
#line 21843 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr672;
		case 32: goto tr671;
		case 46: goto st385;
		case 66: goto st67;
		case 71: goto st389;
		case 78: goto st389;
		case 84: goto st389;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr671;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st364;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st389;
	goto tr377;
tr672:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st365;
tr747:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st365;
st365:
	if ( ++p == pe )
		goto _test_eof365;
case 365:
#line 21994 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr581;
		case 66: goto tr306;
		case 71: goto tr608;
		case 78: goto tr608;
		case 84: goto tr608;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr702;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr608;
	goto tr416;
tr702:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st366;
st366:
	if ( ++p == pe )
		goto _test_eof366;
case 366:
#line 22063 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr704;
		case 32: goto tr703;
		case 46: goto st383;
		case 66: goto st53;
		case 71: goto st317;
		case 78: goto st317;
		case 84: goto st317;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr703;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st366;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st317;
	goto tr416;
tr703:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st367;
tr715:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st367;
st367:
	if ( ++p == pe )
		goto _test_eof367;
case 367:
#line 22170 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr549;
		case 66: goto tr105;
		case 71: goto tr515;
		case 78: goto tr515;
		case 84: goto tr515;
		case 95: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr706;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr105;
		} else if ( (*p) >= 68 )
			goto tr105;
	} else
		goto tr515;
	goto tr421;
tr706:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st368;
st368:
	if ( ++p == pe )
		goto _test_eof368;
case 368:
#line 22218 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr708;
		case 32: goto tr707;
		case 46: goto st277;
		case 66: goto st65;
		case 71: goto st270;
		case 78: goto st270;
		case 84: goto st270;
		case 95: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st368;
		} else if ( (*p) >= 9 )
			goto tr707;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st65;
		} else if ( (*p) >= 68 )
			goto st65;
	} else
		goto st270;
	goto tr421;
tr614:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st369;
tr708:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st369;
st369:
	if ( ++p == pe )
		goto _test_eof369;
case 369:
#line 22345 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr633;
		case 66: goto tr335;
		case 71: goto tr634;
		case 78: goto tr634;
		case 84: goto tr634;
		case 95: goto tr335;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr710;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr335;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr335;
		} else
			goto tr171;
	} else
		goto tr634;
	goto tr426;
tr710:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st370;
st370:
	if ( ++p == pe )
		goto _test_eof370;
case 370:
#line 22410 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr712;
		case 32: goto tr711;
		case 46: goto st387;
		case 66: goto st67;
		case 71: goto st389;
		case 78: goto st389;
		case 84: goto st389;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr711;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st370;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st389;
	goto tr426;
tr636:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st371;
tr712:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st371;
st371:
	if ( ++p == pe )
		goto _test_eof371;
case 371:
#line 22537 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr581;
		case 66: goto tr306;
		case 71: goto tr608;
		case 78: goto tr608;
		case 84: goto tr608;
		case 95: goto tr306;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr714;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr306;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr306;
		} else
			goto tr171;
	} else
		goto tr608;
	goto tr431;
tr714:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st372;
st372:
	if ( ++p == pe )
		goto _test_eof372;
case 372:
#line 22602 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr716;
		case 32: goto tr715;
		case 46: goto st385;
		case 66: goto st53;
		case 71: goto st317;
		case 78: goto st317;
		case 84: goto st317;
		case 95: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr715;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st372;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st53;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st53;
		} else
			goto st13;
	} else
		goto st317;
	goto tr431;
tr704:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st373;
tr716:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st373;
st373:
	if ( ++p == pe )
		goto _test_eof373;
case 373:
#line 22749 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr677;
		case 66: goto tr258;
		case 71: goto tr565;
		case 78: goto tr565;
		case 84: goto tr565;
		case 95: goto tr258;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr718;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr258;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr258;
		} else
			goto tr171;
	} else
		goto tr565;
	goto tr436;
tr718:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st374;
st374:
	if ( ++p == pe )
		goto _test_eof374;
case 374:
#line 22818 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr720;
		case 32: goto tr719;
		case 46: goto st383;
		case 66: goto st50;
		case 71: goto st296;
		case 78: goto st296;
		case 84: goto st296;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr719;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st374;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st50;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st50;
		} else
			goto st13;
	} else
		goto st296;
	goto tr436;
tr693:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st375;
tr720:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st375;
st375:
	if ( ++p == pe )
		goto _test_eof375;
case 375:
#line 22965 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr633;
		case 66: goto tr362;
		case 71: goto tr657;
		case 78: goto tr657;
		case 84: goto tr657;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr722;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr657;
	goto tr441;
tr722:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st376;
st376:
	if ( ++p == pe )
		goto _test_eof376;
case 376:
#line 23034 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr724;
		case 32: goto tr723;
		case 46: goto st387;
		case 66: goto st44;
		case 71: goto st341;
		case 78: goto st341;
		case 84: goto st341;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr723;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st376;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st341;
	goto tr441;
tr659:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st377;
tr724:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st377;
st377:
	if ( ++p == pe )
		goto _test_eof377;
case 377:
#line 23169 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr669;
		case 66: goto tr362;
		case 71: goto tr657;
		case 78: goto tr657;
		case 84: goto tr657;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr726;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr657;
	goto tr446;
tr726:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st378;
st378:
	if ( ++p == pe )
		goto _test_eof378;
case 378:
#line 23238 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr728;
		case 32: goto tr727;
		case 46: goto st385;
		case 66: goto st44;
		case 71: goto st341;
		case 78: goto st341;
		case 84: goto st341;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr727;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st378;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st341;
	goto tr446;
tr727:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st379;
tr735:
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st379;
st379:
	if ( ++p == pe )
		goto _test_eof379;
case 379:
#line 23357 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr573;
		case 66: goto tr85;
		case 71: goto tr569;
		case 78: goto tr569;
		case 84: goto tr569;
		case 95: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto tr730;
		} else if ( (*p) >= 9 )
			goto st97;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr85;
		} else if ( (*p) >= 68 )
			goto tr85;
	} else
		goto tr569;
	goto tr451;
tr730:
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st380;
st380:
	if ( ++p == pe )
		goto _test_eof380;
case 380:
#line 23409 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr732;
		case 32: goto tr731;
		case 46: goto st277;
		case 66: goto st37;
		case 71: goto st298;
		case 78: goto st298;
		case 84: goto st298;
		case 95: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st380;
		} else if ( (*p) >= 9 )
			goto tr731;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st37;
		} else if ( (*p) >= 68 )
			goto st37;
	} else
		goto st298;
	goto tr451;
tr728:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st381;
tr736:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
#line 47 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
#line 59 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
#line 71 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
#line 85 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
#line 99 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
#line 119 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st381;
st381:
	if ( ++p == pe )
		goto _test_eof381;
case 381:
#line 23556 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto tr669;
		case 66: goto tr362;
		case 71: goto tr657;
		case 78: goto tr657;
		case 84: goto tr657;
		case 95: goto tr362;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto st97;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr171;
			} else if ( (*p) >= 48 )
				goto tr734;
		} else
			goto tr171;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto tr362;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr171;
			} else if ( (*p) >= 97 )
				goto tr362;
		} else
			goto tr171;
	} else
		goto tr657;
	goto tr456;
tr734:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 55 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st382;
st382:
	if ( ++p == pe )
		goto _test_eof382;
case 382:
#line 23629 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr736;
		case 32: goto tr735;
		case 46: goto st383;
		case 66: goto st44;
		case 71: goto st341;
		case 78: goto st341;
		case 84: goto st341;
		case 95: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr735;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st382;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st44;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st44;
		} else
			goto st13;
	} else
		goto st341;
	goto tr456;
st383:
	if ( ++p == pe )
		goto _test_eof383;
case 383:
	switch( (*p) ) {
		case 10: goto tr482;
		case 32: goto tr481;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr481;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st384;
	} else
		goto st13;
	goto tr221;
st384:
	if ( ++p == pe )
		goto _test_eof384;
case 384:
	switch( (*p) ) {
		case 10: goto tr740;
		case 32: goto tr739;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr739;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st384;
	} else
		goto st13;
	goto tr221;
st385:
	if ( ++p == pe )
		goto _test_eof385;
case 385:
	switch( (*p) ) {
		case 10: goto tr482;
		case 32: goto tr481;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr481;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st386;
	} else
		goto st13;
	goto tr224;
st386:
	if ( ++p == pe )
		goto _test_eof386;
case 386:
	switch( (*p) ) {
		case 10: goto tr536;
		case 32: goto tr535;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr535;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st386;
	} else
		goto st13;
	goto tr224;
st387:
	if ( ++p == pe )
		goto _test_eof387;
case 387:
	switch( (*p) ) {
		case 10: goto tr482;
		case 32: goto tr481;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr481;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st388;
	} else
		goto st13;
	goto tr195;
st388:
	if ( ++p == pe )
		goto _test_eof388;
case 388:
	switch( (*p) ) {
		case 10: goto tr508;
		case 32: goto tr507;
		case 46: goto st252;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr507;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st388;
	} else
		goto st13;
	goto tr195;
tr634:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st389;
st389:
	if ( ++p == pe )
		goto _test_eof389;
case 389:
#line 23846 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr636;
		case 32: goto tr635;
		case 46: goto st252;
		case 66: goto st67;
		case 71: goto st389;
		case 78: goto st389;
		case 84: goto st389;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr635;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st389;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st389;
	goto tr124;
tr542:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st390;
st390:
	if ( ++p == pe )
		goto _test_eof390;
case 390:
#line 23903 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr544;
		case 32: goto tr543;
		case 46: goto st252;
		case 66: goto st77;
		case 71: goto st390;
		case 78: goto st390;
		case 84: goto st390;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr543;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st390;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st390;
	goto tr144;
tr678:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st391;
st391:
	if ( ++p == pe )
		goto _test_eof391;
case 391:
#line 23968 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr744;
		case 32: goto tr743;
		case 46: goto st383;
		case 66: goto st77;
		case 71: goto st390;
		case 78: goto st390;
		case 84: goto st390;
		case 95: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr743;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st391;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st77;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st77;
		} else
			goto st13;
	} else
		goto st390;
	goto tr385;
tr670:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st392;
st392:
	if ( ++p == pe )
		goto _test_eof392;
case 392:
#line 24037 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr747;
		case 32: goto tr746;
		case 46: goto st383;
		case 66: goto st67;
		case 71: goto st389;
		case 78: goto st389;
		case 84: goto st389;
		case 95: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr746;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st392;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st67;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st67;
		} else
			goto st13;
	} else
		goto st389;
	goto tr374;
st393:
	if ( ++p == pe )
		goto _test_eof393;
case 393:
	switch( (*p) ) {
		case 10: goto tr170;
		case 32: goto st97;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st394;
	} else if ( (*p) >= 9 )
		goto st97;
	goto tr15;
st394:
	if ( ++p == pe )
		goto _test_eof394;
case 394:
	switch( (*p) ) {
		case 10: goto tr499;
		case 32: goto tr498;
		case 46: goto st97;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st394;
	} else if ( (*p) >= 9 )
		goto tr498;
	goto tr15;
tr583:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st395;
st395:
	if ( ++p == pe )
		goto _test_eof395;
case 395:
#line 24134 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr597;
		case 32: goto tr596;
		case 46: goto st252;
		case 66: goto st81;
		case 71: goto st395;
		case 78: goto st395;
		case 84: goto st395;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr596;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st395;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st395;
	goto tr150;
tr489:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st396;
st396:
	if ( ++p == pe )
		goto _test_eof396;
case 396:
#line 24182 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr491;
		case 32: goto tr490;
		case 46: goto st97;
		case 66: goto st73;
		case 71: goto st396;
		case 78: goto st396;
		case 84: goto st396;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st396;
		} else if ( (*p) >= 9 )
			goto tr490;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st396;
	goto tr7;
tr582:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st397;
st397:
	if ( ++p == pe )
		goto _test_eof397;
case 397:
#line 24235 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr751;
		case 32: goto tr750;
		case 46: goto st383;
		case 66: goto st81;
		case 71: goto st395;
		case 78: goto st395;
		case 84: goto st395;
		case 95: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr750;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st13;
			} else if ( (*p) >= 48 )
				goto st397;
		} else
			goto st13;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st81;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st13;
			} else if ( (*p) >= 97 )
				goto st81;
		} else
			goto st13;
	} else
		goto st395;
	goto tr273;
tr574:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st398;
st398:
	if ( ++p == pe )
		goto _test_eof398;
case 398:
#line 24291 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr576;
		case 32: goto tr575;
		case 46: goto st393;
		case 66: goto st26;
		case 71: goto st346;
		case 78: goto st346;
		case 84: goto st346;
		case 95: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st398;
		} else if ( (*p) >= 9 )
			goto tr575;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st26;
		} else if ( (*p) >= 68 )
			goto st26;
	} else
		goto st346;
	goto tr68;
tr550:
#line 67 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st399;
st399:
	if ( ++p == pe )
		goto _test_eof399;
case 399:
#line 24331 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr552;
		case 32: goto tr551;
		case 46: goto st393;
		case 66: goto st73;
		case 71: goto st396;
		case 78: goto st396;
		case 84: goto st396;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st399;
		} else if ( (*p) >= 9 )
			goto tr551;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st73;
		} else if ( (*p) >= 68 )
			goto st73;
	} else
		goto st396;
	goto tr106;
tr534:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 43 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
#line 115 "ped.ragel"
	{
        ts = p;
    }
	goto st400;
st400:
	if ( ++p == pe )
		goto _test_eof400;
case 400:
#line 24380 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr740;
		case 32: goto tr739;
		case 46: goto st383;
		case 65: goto st252;
		case 67: goto st252;
		case 71: goto st252;
		case 78: goto st252;
		case 84: goto st252;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr739;
	} else if ( (*p) > 47 ) {
		if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 126 )
				goto st13;
		} else if ( (*p) >= 48 )
			goto st400;
	} else
		goto st13;
	goto tr221;
tr524:
#line 81 "ped.ragel"
	{
        ts = p;
    }
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st401;
st401:
	if ( ++p == pe )
		goto _test_eof401;
case 401:
#line 24417 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr526;
		case 32: goto tr525;
		case 46: goto st393;
		case 66: goto st25;
		case 71: goto st303;
		case 78: goto st303;
		case 84: goto st303;
		case 95: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( (*p) > 13 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st401;
		} else if ( (*p) >= 9 )
			goto tr525;
	} else if ( (*p) > 67 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 68 )
			goto st25;
	} else
		goto st303;
	goto tr51;
tr497:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st402;
st402:
	if ( ++p == pe )
		goto _test_eof402;
case 402:
#line 24453 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr499;
		case 32: goto tr498;
		case 46: goto st393;
		case 65: goto st97;
		case 67: goto st97;
		case 71: goto st97;
		case 78: goto st97;
		case 84: goto st97;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st402;
	} else if ( (*p) >= 9 )
		goto tr498;
	goto tr15;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st89;
	goto tr19;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	if ( (*p) == 32 )
		goto tr21;
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st89;
	} else if ( (*p) >= 9 )
		goto tr21;
	goto tr19;
tr17:
#line 95 "ped.ragel"
	{
        ts = p;
    }
	goto st90;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
#line 24499 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr18;
		case 46: goto st91;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st90;
	} else if ( (*p) >= 9 )
		goto tr18;
	goto tr15;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st92;
	goto tr15;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	if ( (*p) == 32 )
		goto tr18;
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st92;
	} else if ( (*p) >= 9 )
		goto tr18;
	goto tr15;
tr13:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st93;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
#line 24539 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr14;
		case 95: goto st93;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr14;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st93;
		} else if ( (*p) >= 65 )
			goto st93;
	} else
		goto st93;
	goto tr11;
tr9:
#line 67 "ped.ragel"
	{
        ts = p;
    }
	goto st94;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
#line 24566 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr10;
		case 95: goto st94;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr10;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st94;
		} else if ( (*p) >= 65 )
			goto st94;
	} else
		goto st94;
	goto tr7;
	}
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof153: cs = 153; goto _test_eof; 
	_test_eof154: cs = 154; goto _test_eof; 
	_test_eof155: cs = 155; goto _test_eof; 
	_test_eof156: cs = 156; goto _test_eof; 
	_test_eof157: cs = 157; goto _test_eof; 
	_test_eof158: cs = 158; goto _test_eof; 
	_test_eof159: cs = 159; goto _test_eof; 
	_test_eof160: cs = 160; goto _test_eof; 
	_test_eof161: cs = 161; goto _test_eof; 
	_test_eof162: cs = 162; goto _test_eof; 
	_test_eof163: cs = 163; goto _test_eof; 
	_test_eof164: cs = 164; goto _test_eof; 
	_test_eof165: cs = 165; goto _test_eof; 
	_test_eof166: cs = 166; goto _test_eof; 
	_test_eof167: cs = 167; goto _test_eof; 
	_test_eof168: cs = 168; goto _test_eof; 
	_test_eof169: cs = 169; goto _test_eof; 
	_test_eof170: cs = 170; goto _test_eof; 
	_test_eof171: cs = 171; goto _test_eof; 
	_test_eof172: cs = 172; goto _test_eof; 
	_test_eof173: cs = 173; goto _test_eof; 
	_test_eof174: cs = 174; goto _test_eof; 
	_test_eof175: cs = 175; goto _test_eof; 
	_test_eof176: cs = 176; goto _test_eof; 
	_test_eof177: cs = 177; goto _test_eof; 
	_test_eof178: cs = 178; goto _test_eof; 
	_test_eof179: cs = 179; goto _test_eof; 
	_test_eof180: cs = 180; goto _test_eof; 
	_test_eof181: cs = 181; goto _test_eof; 
	_test_eof182: cs = 182; goto _test_eof; 
	_test_eof183: cs = 183; goto _test_eof; 
	_test_eof184: cs = 184; goto _test_eof; 
	_test_eof185: cs = 185; goto _test_eof; 
	_test_eof186: cs = 186; goto _test_eof; 
	_test_eof187: cs = 187; goto _test_eof; 
	_test_eof188: cs = 188; goto _test_eof; 
	_test_eof189: cs = 189; goto _test_eof; 
	_test_eof190: cs = 190; goto _test_eof; 
	_test_eof191: cs = 191; goto _test_eof; 
	_test_eof192: cs = 192; goto _test_eof; 
	_test_eof193: cs = 193; goto _test_eof; 
	_test_eof194: cs = 194; goto _test_eof; 
	_test_eof195: cs = 195; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof196: cs = 196; goto _test_eof; 
	_test_eof197: cs = 197; goto _test_eof; 
	_test_eof198: cs = 198; goto _test_eof; 
	_test_eof199: cs = 199; goto _test_eof; 
	_test_eof200: cs = 200; goto _test_eof; 
	_test_eof201: cs = 201; goto _test_eof; 
	_test_eof202: cs = 202; goto _test_eof; 
	_test_eof203: cs = 203; goto _test_eof; 
	_test_eof204: cs = 204; goto _test_eof; 
	_test_eof205: cs = 205; goto _test_eof; 
	_test_eof206: cs = 206; goto _test_eof; 
	_test_eof207: cs = 207; goto _test_eof; 
	_test_eof208: cs = 208; goto _test_eof; 
	_test_eof209: cs = 209; goto _test_eof; 
	_test_eof210: cs = 210; goto _test_eof; 
	_test_eof211: cs = 211; goto _test_eof; 
	_test_eof212: cs = 212; goto _test_eof; 
	_test_eof213: cs = 213; goto _test_eof; 
	_test_eof214: cs = 214; goto _test_eof; 
	_test_eof215: cs = 215; goto _test_eof; 
	_test_eof216: cs = 216; goto _test_eof; 
	_test_eof217: cs = 217; goto _test_eof; 
	_test_eof218: cs = 218; goto _test_eof; 
	_test_eof219: cs = 219; goto _test_eof; 
	_test_eof220: cs = 220; goto _test_eof; 
	_test_eof221: cs = 221; goto _test_eof; 
	_test_eof222: cs = 222; goto _test_eof; 
	_test_eof223: cs = 223; goto _test_eof; 
	_test_eof224: cs = 224; goto _test_eof; 
	_test_eof225: cs = 225; goto _test_eof; 
	_test_eof226: cs = 226; goto _test_eof; 
	_test_eof227: cs = 227; goto _test_eof; 
	_test_eof228: cs = 228; goto _test_eof; 
	_test_eof229: cs = 229; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof230: cs = 230; goto _test_eof; 
	_test_eof231: cs = 231; goto _test_eof; 
	_test_eof232: cs = 232; goto _test_eof; 
	_test_eof233: cs = 233; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof234: cs = 234; goto _test_eof; 
	_test_eof235: cs = 235; goto _test_eof; 
	_test_eof236: cs = 236; goto _test_eof; 
	_test_eof237: cs = 237; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof238: cs = 238; goto _test_eof; 
	_test_eof239: cs = 239; goto _test_eof; 
	_test_eof240: cs = 240; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof241: cs = 241; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof242: cs = 242; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof243: cs = 243; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof244: cs = 244; goto _test_eof; 
	_test_eof245: cs = 245; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof246: cs = 246; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof247: cs = 247; goto _test_eof; 
	_test_eof248: cs = 248; goto _test_eof; 
	_test_eof249: cs = 249; goto _test_eof; 
	_test_eof250: cs = 250; goto _test_eof; 
	_test_eof251: cs = 251; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof252: cs = 252; goto _test_eof; 
	_test_eof253: cs = 253; goto _test_eof; 
	_test_eof254: cs = 254; goto _test_eof; 
	_test_eof255: cs = 255; goto _test_eof; 
	_test_eof256: cs = 256; goto _test_eof; 
	_test_eof257: cs = 257; goto _test_eof; 
	_test_eof258: cs = 258; goto _test_eof; 
	_test_eof259: cs = 259; goto _test_eof; 
	_test_eof260: cs = 260; goto _test_eof; 
	_test_eof261: cs = 261; goto _test_eof; 
	_test_eof262: cs = 262; goto _test_eof; 
	_test_eof263: cs = 263; goto _test_eof; 
	_test_eof264: cs = 264; goto _test_eof; 
	_test_eof265: cs = 265; goto _test_eof; 
	_test_eof266: cs = 266; goto _test_eof; 
	_test_eof267: cs = 267; goto _test_eof; 
	_test_eof268: cs = 268; goto _test_eof; 
	_test_eof269: cs = 269; goto _test_eof; 
	_test_eof270: cs = 270; goto _test_eof; 
	_test_eof271: cs = 271; goto _test_eof; 
	_test_eof272: cs = 272; goto _test_eof; 
	_test_eof273: cs = 273; goto _test_eof; 
	_test_eof274: cs = 274; goto _test_eof; 
	_test_eof275: cs = 275; goto _test_eof; 
	_test_eof276: cs = 276; goto _test_eof; 
	_test_eof277: cs = 277; goto _test_eof; 
	_test_eof278: cs = 278; goto _test_eof; 
	_test_eof279: cs = 279; goto _test_eof; 
	_test_eof280: cs = 280; goto _test_eof; 
	_test_eof281: cs = 281; goto _test_eof; 
	_test_eof282: cs = 282; goto _test_eof; 
	_test_eof283: cs = 283; goto _test_eof; 
	_test_eof284: cs = 284; goto _test_eof; 
	_test_eof285: cs = 285; goto _test_eof; 
	_test_eof286: cs = 286; goto _test_eof; 
	_test_eof287: cs = 287; goto _test_eof; 
	_test_eof288: cs = 288; goto _test_eof; 
	_test_eof289: cs = 289; goto _test_eof; 
	_test_eof290: cs = 290; goto _test_eof; 
	_test_eof291: cs = 291; goto _test_eof; 
	_test_eof292: cs = 292; goto _test_eof; 
	_test_eof293: cs = 293; goto _test_eof; 
	_test_eof294: cs = 294; goto _test_eof; 
	_test_eof295: cs = 295; goto _test_eof; 
	_test_eof296: cs = 296; goto _test_eof; 
	_test_eof297: cs = 297; goto _test_eof; 
	_test_eof298: cs = 298; goto _test_eof; 
	_test_eof299: cs = 299; goto _test_eof; 
	_test_eof300: cs = 300; goto _test_eof; 
	_test_eof301: cs = 301; goto _test_eof; 
	_test_eof302: cs = 302; goto _test_eof; 
	_test_eof303: cs = 303; goto _test_eof; 
	_test_eof304: cs = 304; goto _test_eof; 
	_test_eof305: cs = 305; goto _test_eof; 
	_test_eof306: cs = 306; goto _test_eof; 
	_test_eof307: cs = 307; goto _test_eof; 
	_test_eof308: cs = 308; goto _test_eof; 
	_test_eof309: cs = 309; goto _test_eof; 
	_test_eof310: cs = 310; goto _test_eof; 
	_test_eof311: cs = 311; goto _test_eof; 
	_test_eof312: cs = 312; goto _test_eof; 
	_test_eof313: cs = 313; goto _test_eof; 
	_test_eof314: cs = 314; goto _test_eof; 
	_test_eof315: cs = 315; goto _test_eof; 
	_test_eof316: cs = 316; goto _test_eof; 
	_test_eof317: cs = 317; goto _test_eof; 
	_test_eof318: cs = 318; goto _test_eof; 
	_test_eof319: cs = 319; goto _test_eof; 
	_test_eof320: cs = 320; goto _test_eof; 
	_test_eof321: cs = 321; goto _test_eof; 
	_test_eof322: cs = 322; goto _test_eof; 
	_test_eof323: cs = 323; goto _test_eof; 
	_test_eof324: cs = 324; goto _test_eof; 
	_test_eof325: cs = 325; goto _test_eof; 
	_test_eof326: cs = 326; goto _test_eof; 
	_test_eof327: cs = 327; goto _test_eof; 
	_test_eof328: cs = 328; goto _test_eof; 
	_test_eof329: cs = 329; goto _test_eof; 
	_test_eof330: cs = 330; goto _test_eof; 
	_test_eof331: cs = 331; goto _test_eof; 
	_test_eof332: cs = 332; goto _test_eof; 
	_test_eof333: cs = 333; goto _test_eof; 
	_test_eof334: cs = 334; goto _test_eof; 
	_test_eof335: cs = 335; goto _test_eof; 
	_test_eof336: cs = 336; goto _test_eof; 
	_test_eof337: cs = 337; goto _test_eof; 
	_test_eof338: cs = 338; goto _test_eof; 
	_test_eof339: cs = 339; goto _test_eof; 
	_test_eof340: cs = 340; goto _test_eof; 
	_test_eof341: cs = 341; goto _test_eof; 
	_test_eof342: cs = 342; goto _test_eof; 
	_test_eof343: cs = 343; goto _test_eof; 
	_test_eof344: cs = 344; goto _test_eof; 
	_test_eof345: cs = 345; goto _test_eof; 
	_test_eof346: cs = 346; goto _test_eof; 
	_test_eof347: cs = 347; goto _test_eof; 
	_test_eof348: cs = 348; goto _test_eof; 
	_test_eof349: cs = 349; goto _test_eof; 
	_test_eof350: cs = 350; goto _test_eof; 
	_test_eof351: cs = 351; goto _test_eof; 
	_test_eof352: cs = 352; goto _test_eof; 
	_test_eof353: cs = 353; goto _test_eof; 
	_test_eof354: cs = 354; goto _test_eof; 
	_test_eof355: cs = 355; goto _test_eof; 
	_test_eof356: cs = 356; goto _test_eof; 
	_test_eof357: cs = 357; goto _test_eof; 
	_test_eof358: cs = 358; goto _test_eof; 
	_test_eof359: cs = 359; goto _test_eof; 
	_test_eof360: cs = 360; goto _test_eof; 
	_test_eof361: cs = 361; goto _test_eof; 
	_test_eof362: cs = 362; goto _test_eof; 
	_test_eof363: cs = 363; goto _test_eof; 
	_test_eof364: cs = 364; goto _test_eof; 
	_test_eof365: cs = 365; goto _test_eof; 
	_test_eof366: cs = 366; goto _test_eof; 
	_test_eof367: cs = 367; goto _test_eof; 
	_test_eof368: cs = 368; goto _test_eof; 
	_test_eof369: cs = 369; goto _test_eof; 
	_test_eof370: cs = 370; goto _test_eof; 
	_test_eof371: cs = 371; goto _test_eof; 
	_test_eof372: cs = 372; goto _test_eof; 
	_test_eof373: cs = 373; goto _test_eof; 
	_test_eof374: cs = 374; goto _test_eof; 
	_test_eof375: cs = 375; goto _test_eof; 
	_test_eof376: cs = 376; goto _test_eof; 
	_test_eof377: cs = 377; goto _test_eof; 
	_test_eof378: cs = 378; goto _test_eof; 
	_test_eof379: cs = 379; goto _test_eof; 
	_test_eof380: cs = 380; goto _test_eof; 
	_test_eof381: cs = 381; goto _test_eof; 
	_test_eof382: cs = 382; goto _test_eof; 
	_test_eof383: cs = 383; goto _test_eof; 
	_test_eof384: cs = 384; goto _test_eof; 
	_test_eof385: cs = 385; goto _test_eof; 
	_test_eof386: cs = 386; goto _test_eof; 
	_test_eof387: cs = 387; goto _test_eof; 
	_test_eof388: cs = 388; goto _test_eof; 
	_test_eof389: cs = 389; goto _test_eof; 
	_test_eof390: cs = 390; goto _test_eof; 
	_test_eof391: cs = 391; goto _test_eof; 
	_test_eof392: cs = 392; goto _test_eof; 
	_test_eof393: cs = 393; goto _test_eof; 
	_test_eof394: cs = 394; goto _test_eof; 
	_test_eof395: cs = 395; goto _test_eof; 
	_test_eof396: cs = 396; goto _test_eof; 
	_test_eof397: cs = 397; goto _test_eof; 
	_test_eof398: cs = 398; goto _test_eof; 
	_test_eof399: cs = 399; goto _test_eof; 
	_test_eof400: cs = 400; goto _test_eof; 
	_test_eof401: cs = 401; goto _test_eof; 
	_test_eof402: cs = 402; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 97: 
	case 98: 
	case 99: 
	case 100: 
	case 101: 
	case 102: 
	case 103: 
	case 104: 
	case 105: 
	case 106: 
	case 107: 
	case 108: 
	case 109: 
	case 110: 
	case 111: 
	case 112: 
	case 113: 
	case 114: 
	case 115: 
	case 116: 
	case 117: 
	case 118: 
	case 119: 
	case 120: 
	case 121: 
	case 122: 
	case 123: 
	case 124: 
	case 125: 
	case 126: 
	case 127: 
	case 128: 
	case 129: 
	case 130: 
	case 131: 
	case 132: 
	case 133: 
	case 134: 
	case 135: 
	case 136: 
	case 137: 
	case 138: 
	case 139: 
	case 140: 
	case 141: 
	case 142: 
	case 143: 
	case 144: 
	case 145: 
	case 146: 
	case 147: 
	case 148: 
	case 149: 
	case 150: 
	case 151: 
	case 152: 
	case 153: 
	case 154: 
	case 155: 
	case 156: 
	case 157: 
	case 158: 
	case 159: 
	case 160: 
	case 161: 
	case 162: 
	case 163: 
	case 164: 
	case 165: 
	case 166: 
	case 167: 
	case 168: 
	case 169: 
	case 170: 
	case 171: 
	case 172: 
	case 173: 
	case 174: 
	case 175: 
	case 176: 
	case 177: 
	case 178: 
	case 179: 
	case 180: 
	case 181: 
	case 182: 
	case 183: 
	case 184: 
	case 185: 
	case 186: 
	case 187: 
	case 188: 
	case 189: 
	case 190: 
	case 191: 
	case 192: 
	case 193: 
	case 194: 
	case 195: 
	case 196: 
	case 197: 
	case 198: 
	case 199: 
	case 200: 
	case 201: 
	case 202: 
	case 203: 
	case 204: 
	case 205: 
	case 206: 
	case 207: 
	case 208: 
	case 209: 
	case 210: 
	case 211: 
	case 212: 
	case 213: 
	case 214: 
	case 215: 
	case 216: 
	case 217: 
	case 218: 
	case 219: 
	case 220: 
	case 221: 
	case 222: 
	case 223: 
	case 224: 
	case 225: 
	case 226: 
	case 227: 
	case 228: 
	case 229: 
	case 230: 
	case 231: 
	case 232: 
	case 233: 
	case 234: 
	case 235: 
	case 236: 
	case 237: 
	case 238: 
	case 239: 
	case 240: 
	case 241: 
	case 242: 
	case 243: 
	case 244: 
	case 245: 
	case 246: 
	case 247: 
	case 248: 
	case 249: 
	case 250: 
	case 251: 
	case 252: 
	case 253: 
	case 254: 
	case 255: 
	case 256: 
	case 257: 
	case 258: 
	case 259: 
	case 260: 
	case 261: 
	case 262: 
	case 263: 
	case 264: 
	case 265: 
	case 266: 
	case 267: 
	case 268: 
	case 269: 
	case 270: 
	case 271: 
	case 272: 
	case 273: 
	case 274: 
	case 275: 
	case 276: 
	case 277: 
	case 278: 
	case 279: 
	case 280: 
	case 281: 
	case 282: 
	case 283: 
	case 284: 
	case 285: 
	case 286: 
	case 287: 
	case 288: 
	case 289: 
	case 290: 
	case 291: 
	case 292: 
	case 293: 
	case 294: 
	case 295: 
	case 296: 
	case 297: 
	case 298: 
	case 299: 
	case 300: 
	case 301: 
	case 302: 
	case 303: 
	case 304: 
	case 305: 
	case 306: 
	case 307: 
	case 308: 
	case 309: 
	case 310: 
	case 311: 
	case 312: 
	case 313: 
	case 314: 
	case 315: 
	case 316: 
	case 317: 
	case 318: 
	case 319: 
	case 320: 
	case 321: 
	case 322: 
	case 323: 
	case 324: 
	case 325: 
	case 326: 
	case 327: 
	case 328: 
	case 329: 
	case 330: 
	case 331: 
	case 332: 
	case 333: 
	case 334: 
	case 335: 
	case 336: 
	case 337: 
	case 338: 
	case 339: 
	case 340: 
	case 341: 
	case 342: 
	case 343: 
	case 344: 
	case 345: 
	case 346: 
	case 347: 
	case 348: 
	case 349: 
	case 350: 
	case 351: 
	case 352: 
	case 353: 
	case 354: 
	case 355: 
	case 356: 
	case 357: 
	case 358: 
	case 359: 
	case 360: 
	case 361: 
	case 362: 
	case 363: 
	case 364: 
	case 365: 
	case 366: 
	case 367: 
	case 368: 
	case 369: 
	case 370: 
	case 371: 
	case 372: 
	case 373: 
	case 374: 
	case 375: 
	case 376: 
	case 377: 
	case 378: 
	case 379: 
	case 380: 
	case 381: 
	case 382: 
	case 383: 
	case 384: 
	case 385: 
	case 386: 
	case 387: 
	case 388: 
	case 389: 
	case 390: 
	case 391: 
	case 392: 
	case 393: 
	case 394: 
	case 395: 
	case 396: 
	case 397: 
	case 398: 
	case 399: 
	case 400: 
	case 401: 
	case 402: 
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        add_record_to_ped_batch(current_record, current_batch);

        num_records++;
    }
	break;
	case 1: 
	case 13: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	break;
	case 2: 
	case 3: 
	case 14: 
	case 15: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	break;
	case 4: 
	case 5: 
	case 16: 
	case 17: 
	case 73: 
	case 94: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 18: 
	case 19: 
	case 25: 
	case 86: 
	case 93: 
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 20: 
	case 21: 
	case 35: 
	case 36: 
	case 58: 
	case 87: 
	case 90: 
	case 91: 
	case 92: 
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 10: 
	case 11: 
	case 22: 
	case 23: 
	case 31: 
	case 59: 
	case 60: 
	case 88: 
	case 89: 
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 79: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	break;
	case 77: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 81: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 65: 
	case 80: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 74: 
	case 78: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 82: 
	case 83: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 26: 
	case 52: 
	case 66: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 55: 
	case 75: 
	case 76: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 84: 
	case 85: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 27: 
	case 28: 
	case 34: 
	case 64: 
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 56: 
	case 57: 
	case 72: 
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 29: 
	case 30: 
	case 32: 
	case 33: 
	case 42: 
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 50: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 53: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 67: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 37: 
	case 51: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 54: 
	case 61: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 68: 
	case 69: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 38: 
	case 39: 
	case 43: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 62: 
	case 63: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 70: 
	case 71: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 40: 
	case 41: 
	case 49: 
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 44: 
#line 51 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 45: 
	case 46: 
#line 63 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 47: 
	case 48: 
#line 77 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
#line 91 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
#line 111 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
#line 25724 "ped_reader.c"
	}
	}

	_out: {}
	}

#line 182 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 25743 "ped_reader.c"
95
#line 192 "ped.ragel"
 ) 
    {
        LOG_ERROR("The file was not successfully read\n");
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 25751 "ped_reader.c"
95
#line 196 "ped.ragel"
);
    } 

    LOG_INFO_F("Records read = %zu\n", num_records);

    return cs < 
#line 25760 "ped_reader.c"
95
#line 201 "ped.ragel"
;
}
