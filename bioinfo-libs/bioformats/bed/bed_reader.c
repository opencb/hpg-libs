
#line 1 "bed.ragel"
#include "bed_reader.h"

static size_t lines = 1;
static size_t num_records = 0;
static size_t num_batches = 0;

static bed_record_t *current_record;
static bed_header_entry_t *current_header_entry;
static bed_batch_t *current_batch;


#line 15 "bed_reader.c"
static const int bed_start = 62;
static const int bed_first_final = 62;
static const int bed_error = 0;

static const int bed_en_main = 62;


#line 273 "bed.ragel"



int bed_ragel_read(list_t *batches_list, size_t batch_size, bed_file_t *file) {
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_header_entry = bed_header_entry_new();
    current_batch = bed_batch_new(batch_size);

    
#line 40 "bed_reader.c"
	{
	cs = bed_start;
	}

#line 45 "bed_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 62:
	switch( (*p) ) {
		case 10: goto tr85;
		case 99: goto tr86;
		case 115: goto tr87;
	}
	if ( 0 <= (*p) )
		goto tr84;
	goto tr3;
tr0:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	goto st0;
tr3:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'chrom' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr8:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 85 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr13:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 99 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr18:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr22:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr26:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr29:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr33:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr37:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr41:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr45:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr49:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
#line 186 "bed_reader.c"
st0:
cs = 0;
	goto _out;
tr84:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 201 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( 0 <= (*p) )
		goto st1;
	goto tr0;
tr2:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st63;
tr88:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
#line 240 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr88;
		case 99: goto tr86;
		case 115: goto tr87;
	}
	if ( 0 <= (*p) )
		goto tr84;
	goto tr3;
tr86:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st2;
tr94:
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 278 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 104: goto st3;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 10: goto tr2;
		case 114: goto st4;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 10: goto tr2;
		case 95: goto st5;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st5;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st5;
	} else
		goto st5;
	goto tr3;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 9: goto tr7;
		case 10: goto tr2;
		case 95: goto st5;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st5;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st5;
	} else
		goto st5;
	goto tr3;
tr7:
#line 67 "bed.ragel"
	{
        set_bed_record_chromosome(ts, p-ts, current_record);
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 366 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr9;
	goto tr8;
tr9:
#line 75 "bed.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 388 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 10: goto tr2;
		case 46: goto st60;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st7;
	goto tr8;
tr10:
#line 79 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_start(atol(field), current_record);
        free(field);
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 415 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr14;
	goto tr13;
tr14:
#line 89 "bed.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 437 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr15;
		case 10: goto tr2;
		case 46: goto st58;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st9;
	goto tr13;
tr15:
#line 93 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_end(atol(field), current_record);
        free(field);
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 464 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr19;
		case 95: goto tr20;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto tr20;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto tr20;
	} else
		goto tr20;
	goto tr18;
tr19:
#line 103 "bed.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 501 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr21;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr18;
tr21:
#line 107 "bed.ragel"
	{
        set_bed_record_name(ts, p-ts, current_record);
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
#line 519 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr23;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr24;
	goto tr22;
tr23:
#line 115 "bed.ragel"
	{
        ts = p;
    }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 543 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr25;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr22;
tr25:
#line 119 "bed.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_bed_record_score(score, current_record);
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 567 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 43: goto tr27;
	}
	if ( (*p) < 45 ) {
		if ( 0 <= (*p) && (*p) <= 44 )
			goto st1;
	} else if ( (*p) > 46 ) {
		if ( 47 <= (*p) )
			goto st1;
	} else
		goto tr27;
	goto tr26;
tr27:
#line 133 "bed.ragel"
	{
        ts = p;
    }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 591 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr28;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr26;
tr28:
#line 137 "bed.ragel"
	{
        set_bed_record_strand(*ts, current_record);
    }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 609 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr30;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr31;
	goto tr29;
tr30:
#line 145 "bed.ragel"
	{
        ts = p;
    }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
#line 633 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr32;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr29;
tr32:
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 653 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr34;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr35;
	goto tr33;
tr34:
#line 159 "bed.ragel"
	{
        ts = p;
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 677 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr36;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr33;
tr36:
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 697 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr38;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr39;
	goto tr37;
tr38:
#line 173 "bed.ragel"
	{
        ts = p;
    }
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
#line 721 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr40;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr37;
tr40:
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
#line 739 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr42;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr43;
	goto tr41;
tr42:
#line 185 "bed.ragel"
	{
        ts = p;
    }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 763 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr44;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr41;
tr44:
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
#line 783 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr46;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr47;
	goto tr45;
tr46:
#line 199 "bed.ragel"
	{
        ts = p;
    }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 807 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr48;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr45;
tr48:
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 825 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 46: goto tr50;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr51;
	goto tr49;
tr50:
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 849 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr89;
	if ( 0 <= (*p) )
		goto st1;
	goto tr49;
tr91:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st65;
tr89:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 215 "bed.ragel"
	{
        set_bed_record_blockstarts(ts, p-ts, current_record);
    }
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 908 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr88;
		case 35: goto tr90;
		case 99: goto tr86;
		case 115: goto tr87;
	}
	if ( 0 <= (*p) )
		goto tr84;
	goto tr3;
tr90:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 929 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 127: goto st1;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st66;
	} else if ( (*p) >= 0 )
		goto st1;
	goto tr0;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 10: goto tr91;
		case 127: goto st1;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st66;
	} else if ( (*p) >= 0 )
		goto st1;
	goto tr0;
tr87:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st28;
tr95:
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st28;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
#line 983 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto st29;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 10: goto tr2;
		case 97: goto st30;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 10: goto tr2;
		case 102: goto st31;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 10: goto tr2;
		case 102: goto st32;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 10: goto tr2;
		case 111: goto st33;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 10: goto tr2;
		case 108: goto st34;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 10: goto tr2;
		case 100: goto st4;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
tr51:
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 1067 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr89;
		case 44: goto st35;
		case 46: goto st36;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st67;
	goto tr49;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st67;
	goto tr49;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st68;
	goto tr49;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 10: goto tr89;
		case 44: goto st35;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st68;
	goto tr49;
tr47:
#line 199 "bed.ragel"
	{
        ts = p;
    }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 1139 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr48;
		case 10: goto tr2;
		case 44: goto st38;
		case 46: goto st39;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st37;
	goto tr45;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st37;
	goto tr45;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st40;
	goto tr45;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 9: goto tr48;
		case 10: goto tr2;
		case 44: goto st38;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st40;
	goto tr45;
tr43:
#line 185 "bed.ragel"
	{
        ts = p;
    }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
#line 1213 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr44;
		case 10: goto tr2;
		case 46: goto st42;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st41;
	goto tr41;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st43;
	goto tr41;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 9: goto tr44;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st43;
	goto tr41;
tr39:
#line 173 "bed.ragel"
	{
        ts = p;
    }
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
#line 1270 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr40;
		case 10: goto tr2;
		case 44: goto st45;
		case 46: goto st46;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st44;
	goto tr37;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st44;
	goto tr37;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st47;
	goto tr37;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 9: goto tr40;
		case 10: goto tr2;
		case 44: goto st45;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st47;
	goto tr37;
tr35:
#line 159 "bed.ragel"
	{
        ts = p;
    }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 1344 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr36;
		case 10: goto tr2;
		case 46: goto st49;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st48;
	goto tr33;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st50;
	goto tr33;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 9: goto tr36;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st50;
	goto tr33;
tr31:
#line 145 "bed.ragel"
	{
        ts = p;
    }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
#line 1401 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr32;
		case 10: goto tr2;
		case 46: goto st52;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st51;
	goto tr29;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st53;
	goto tr29;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 9: goto tr32;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st53;
	goto tr29;
tr24:
#line 115 "bed.ragel"
	{
        ts = p;
    }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 1458 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr25;
		case 10: goto tr2;
		case 46: goto st55;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st54;
	goto tr22;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st56;
	goto tr22;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 9: goto tr25;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st56;
	goto tr22;
tr20:
#line 103 "bed.ragel"
	{
        ts = p;
    }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 1515 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr21;
		case 10: goto tr2;
		case 95: goto st57;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st57;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st57;
	} else
		goto st57;
	goto tr18;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st59;
	goto tr13;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 9: goto tr15;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st59;
	goto tr13;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st61;
	goto tr8;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 9: goto tr10;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st61;
	goto tr8;
tr85:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
#line 1617 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto tr94;
		case 115: goto tr95;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
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
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: 
	case 27: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	break;
	case 2: 
	case 3: 
	case 4: 
	case 5: 
	case 28: 
	case 29: 
	case 30: 
	case 31: 
	case 32: 
	case 33: 
	case 34: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'chrom' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 60: 
	case 61: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 85 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 58: 
	case 59: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 99 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 10: 
	case 11: 
	case 57: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 12: 
	case 13: 
	case 54: 
	case 55: 
	case 56: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 14: 
	case 15: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 16: 
	case 17: 
	case 51: 
	case 52: 
	case 53: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 18: 
	case 19: 
	case 48: 
	case 49: 
	case 50: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 20: 
	case 21: 
	case 44: 
	case 45: 
	case 46: 
	case 47: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 22: 
	case 23: 
	case 41: 
	case 42: 
	case 43: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 24: 
	case 25: 
	case 37: 
	case 38: 
	case 39: 
	case 40: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 26: 
	case 35: 
	case 36: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 64: 
	case 67: 
	case 68: 
#line 215 "bed.ragel"
	{
        set_bed_record_blockstarts(ts, p-ts, current_record);
    }
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	break;
#line 1903 "bed_reader.c"
	}
	}

	_out: {}
	}

#line 291 "bed.ragel"
 

    // Insert the last batch
    if (!bed_batch_is_empty(current_batch)) {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->records->size);
    }

    if ( cs < 
#line 1921 "bed_reader.c"
62
#line 300 "bed.ragel"
 ) {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 1927 "bed_reader.c"
62
#line 302 "bed.ragel"
);
    } 

    LOG_INFO_F("BED records read = %zu\n", num_batches * batch_size + num_records);

    // Free current_xxx pointers if not needed in another module
    //bed_header_entry_free(current_header_entry);

    return cs < 
#line 1939 "bed_reader.c"
62
#line 310 "bed.ragel"
;
}
