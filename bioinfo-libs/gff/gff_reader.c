
#line 1 "gff.ragel"
#include "gff_reader.h"

static size_t records = 0;

static gff_record_t *current_record;
static gff_header_entry_t *current_header_entry;
static gff_batch_t *current_batch;

static enum GFF_Field current_field = SEQUENCE;


#line 15 "gff_reader.c"
static const int gff_start = 1;
static const int gff_first_final = 29;
static const int gff_error = 0;

static const int gff_en_header_scan = 32;
static const int gff_en_record_scan = 34;
static const int gff_en_main = 1;


#line 142 "gff.ragel"



static char* get_token(char *ts, char *te)
{
    char *field = (char*) malloc ((te-ts+1) * sizeof(char));
    strncpy(field, ts, (te-ts));
    field[te-ts] = '\0';
    return field;
}

static void set_field(char* ts, char *te)
{
    char *field = get_token(ts, te);
    float float_val = -1.0f;

    switch (current_field)
    {
        case SEQUENCE:
            current_record = create_gff_record();
            set_gff_record_sequence(field, current_record);
        break;
        case SOURCE:
            set_gff_record_source(field, current_record);
        break;
        case FEATURE:
            set_gff_record_feature(field, current_record);
        break;
        case START:
            set_gff_record_start(atol(field), current_record);
            free(field);    // Not set as gff_record_t variable -> not freed later
        break;
        case END:
            set_gff_record_end(atol(field), current_record);
            free(field);    // Not set as gff_record_t variable -> not freed later
        break;
        case SCORE:
            if (strncmp(".", field, 1) != 0)
            {
                float_val = atof(field);
            }
            set_gff_record_score(float_val, current_record);
            free(field);	// Not set as gff_record_t variable -> not freed later
        break;
        case STRAND:
            if (field != NULL && strlen(field) > 0) {
                set_gff_record_strand(field[0], current_record);
            }
        break;
        case FRAME:
            if (strncmp(".", field, 1) != 0)
            {
                float_val = atof(field);
            }
            set_gff_record_frame(float_val, current_record);
            free(field);    // Not set as gff_record_t variable -> not freed later
        break;
        case ATTRIBUTE:
            set_gff_record_attribute(field, current_record);
        break;
    }

    current_field++;
}

int gff_ragel_read(list_t *batches_list, size_t batch_size, gff_file_t *file)
{
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_header_entry = create_gff_header_entry();
    current_batch = gff_batch_new(batch_size);

    
#line 105 "gff_reader.c"
	{
	cs = gff_start;
	top = 0;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 114 "gff_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	goto _resume;

_again:
	switch ( cs ) {
		case 1: goto st1;
		case 0: goto st0;
		case 2: goto st2;
		case 3: goto st3;
		case 4: goto st4;
		case 5: goto st5;
		case 6: goto st6;
		case 7: goto st7;
		case 8: goto st8;
		case 9: goto st9;
		case 10: goto st10;
		case 11: goto st11;
		case 12: goto st12;
		case 13: goto st13;
		case 14: goto st14;
		case 15: goto st15;
		case 16: goto st16;
		case 17: goto st17;
		case 29: goto st29;
		case 30: goto st30;
		case 18: goto st18;
		case 31: goto st31;
		case 19: goto st19;
		case 20: goto st20;
		case 21: goto st21;
		case 22: goto st22;
		case 23: goto st23;
		case 24: goto st24;
		case 25: goto st25;
		case 32: goto st32;
		case 33: goto st33;
		case 34: goto st34;
		case 35: goto st35;
		case 36: goto st36;
		case 26: goto st26;
		case 37: goto st37;
		case 27: goto st27;
		case 38: goto st38;
		case 28: goto st28;
		case 39: goto st39;
		case 40: goto st40;
	default: break;
	}

	if ( ++p == pe )
		goto _test_eof;
_resume:
	switch ( cs )
	{
st1:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 177 "gff_reader.c"
	if ( (*p) == 35 )
		goto st2;
	goto st0;
st0:
cs = 0;
	goto _out;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( (*p) == 35 )
		goto tr2;
	goto st0;
tr2:
#line 27 "gff.ragel"
	{ {stack[top++] = 3; goto st32;} }
	goto st3;
st3:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof3;
case 3:
#line 201 "gff_reader.c"
	switch( (*p) ) {
		case 35: goto st2;
		case 95: goto tr3;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr3;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr3;
	} else
		goto tr3;
	goto st0;
tr3:
#line 140 "gff.ragel"
	{ p--; {stack[top++] = 4; goto st34;} }
	goto st4;
st4:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 225 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto st5;
		case 95: goto st4;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st4;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st4;
	} else
		goto st4;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 95 )
		goto st6;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st6;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st6;
	} else
		goto st6;
	goto st0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	switch( (*p) ) {
		case 9: goto st7;
		case 95: goto st6;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st6;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st6;
	} else
		goto st6;
	goto st0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 95 )
		goto st8;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st8;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st8;
	} else
		goto st8;
	goto st0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	switch( (*p) ) {
		case 9: goto st9;
		case 95: goto st8;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st8;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st8;
	} else
		goto st8;
	goto st0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st10;
	goto st0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	switch( (*p) ) {
		case 9: goto st11;
		case 46: goto st24;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st10;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 9: goto st13;
		case 46: goto st22;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st12;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 46 )
		goto st14;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 9 )
		goto st15;
	goto st0;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 43 )
		goto st16;
	if ( 45 <= (*p) && (*p) <= 46 )
		goto st16;
	goto st0;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( (*p) == 9 )
		goto st17;
	goto st0;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 95 )
		goto st29;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st29;
	} else
		goto st29;
	goto st0;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 10: goto st30;
		case 59: goto st17;
		case 61: goto st18;
		case 95: goto st29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st29;
	} else
		goto st29;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 95 )
		goto st4;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st4;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st4;
	} else
		goto st4;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 10 )
		goto st30;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st31;
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 9: goto st15;
		case 46: goto st20;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 9 )
		goto st15;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st21;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 9 )
		goto st13;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st23;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st25;
	goto st0;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 9 )
		goto st11;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st25;
	goto st0;
tr35:
#line 36 "gff.ragel"
	{te = p+1;{
            add_gff_header_entry(current_header_entry, file);
            current_header_entry = create_gff_header_entry();
            LOG_DEBUG("\n");
            {cs = stack[--top];goto _again;}
        }}
	goto st32;
tr36:
#line 31 "gff.ragel"
	{te = p;p--;{
            char *text = get_token(ts, te);
            set_gff_header_entry_text(text, current_header_entry);
        }}
	goto st32;
st32:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof32;
case 32:
#line 1 "NONE"
	{ts = p;}
#line 519 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr35;
	if ( (*p) > 13 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st33;
	} else if ( (*p) >= 9 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 9 )
		goto st33;
	if ( (*p) > 13 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st33;
	} else if ( (*p) >= 11 )
		goto st33;
	goto tr36;
tr27:
#line 75 "gff.ragel"
	{{p = ((te))-1;}{
            set_field(ts, te);
        }}
	goto st34;
tr29:
#line 1 "NONE"
	{	switch( act ) {
	case 4:
	{{p = ((te))-1;}
            set_field(ts, te);
        }
	break;
	case 12:
	{{p = ((te))-1;}
            set_field(ts, te);
        }
	break;
	}
	}
	goto st34;
tr37:
#line 132 "gff.ragel"
	{te = p+1;}
	goto st34;
tr38:
#line 112 "gff.ragel"
	{te = p+1;{
            // If batch is full, add to the list of batches and create a new, empty one
            if (gff_batch_is_full(current_batch))
            {
                list_item_t *item = list_item_new(records, 1, current_batch); 
                list_insert_item(item, batches_list);
                LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
                current_batch = gff_batch_new(batch_size);
            }
            
            // If not a blank line, add current record to current batch
            if (current_field != SEQUENCE) {
                add_record_to_gff_batch(current_record, current_batch);
            }
            
            current_field = SEQUENCE;
            records++;
            LOG_DEBUG("\n");
        }}
	goto st34;
tr40:
#line 99 "gff.ragel"
	{te = p+1;{
            set_field(ts, te);
        }}
	goto st34;
tr41:
#line 95 "gff.ragel"
	{te = p+1;{
            set_field(ts, te);
        }}
	goto st34;
tr44:
#line 71 "gff.ragel"
	{te = p;p--;{
            LOG_DEBUG("Comment found, nothing to do.");
        }}
	goto st34;
tr45:
#line 75 "gff.ragel"
	{te = p;p--;{
            set_field(ts, te);
        }}
	goto st34;
tr49:
#line 87 "gff.ragel"
	{te = p;p--;{
            set_field(ts, te);
        }}
	goto st34;
tr50:
#line 107 "gff.ragel"
	{te = p;p--;{
            set_field(ts, te);
        }}
	goto st34;
st34:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof34;
case 34:
#line 1 "NONE"
	{ts = p;}
#line 632 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr38;
		case 32: goto tr37;
		case 35: goto st35;
		case 43: goto tr40;
		case 45: goto tr40;
		case 46: goto tr41;
		case 95: goto tr43;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto tr37;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto tr43;
		} else if ( (*p) >= 65 )
			goto tr43;
	} else
		goto tr42;
	goto st0;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st35;
	goto tr44;
tr42:
#line 1 "NONE"
	{te = p+1;}
#line 75 "gff.ragel"
	{act = 4;}
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
#line 671 "gff_reader.c"
	switch( (*p) ) {
		case 46: goto st26;
		case 59: goto st27;
		case 61: goto st28;
		case 95: goto tr43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr42;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr43;
	} else
		goto tr43;
	goto tr45;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st37;
	goto tr27;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st37;
	goto tr49;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 95 )
		goto tr30;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr30;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr30;
	} else
		goto tr30;
	goto tr29;
tr30:
#line 1 "NONE"
	{te = p+1;}
#line 107 "gff.ragel"
	{act = 12;}
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 726 "gff_reader.c"
	switch( (*p) ) {
		case 59: goto st27;
		case 61: goto st28;
		case 95: goto tr30;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr30;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr30;
	} else
		goto tr30;
	goto tr50;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st39;
	goto tr29;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st39;
	goto tr50;
tr43:
#line 1 "NONE"
	{te = p+1;}
#line 75 "gff.ragel"
	{act = 4;}
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
#line 765 "gff_reader.c"
	switch( (*p) ) {
		case 59: goto st27;
		case 61: goto st28;
		case 95: goto tr43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr43;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr43;
	} else
		goto tr43;
	goto tr45;
	}
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
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 33: goto tr36;
	case 35: goto tr44;
	case 36: goto tr45;
	case 26: goto tr27;
	case 37: goto tr49;
	case 27: goto tr29;
	case 38: goto tr50;
	case 28: goto tr29;
	case 39: goto tr50;
	case 40: goto tr45;
	}
	}

	_out: {}
	}

#line 223 "gff.ragel"
 

    // Insert the last batch
    if (!gff_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 854 "gff_reader.c"
29
#line 233 "gff.ragel"
 ) 
    {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 861 "gff_reader.c"
29
#line 236 "gff.ragel"
);
    } 

    LOG_INFO_F("Records read = %zu\n", records);

    // Free current_xxx pointers if not needed in another module
    gff_header_entry_free(current_header_entry);

    return cs < 
#line 873 "gff_reader.c"
29
#line 244 "gff.ragel"
;
}
