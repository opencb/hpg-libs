
#line 1 "ped.ragel"
#include "ped_reader.h"

static size_t lines = 0;
static size_t records = 0;

static ped_record_t *current_record;
static ped_batch_t *current_batch;

static int record_error = 0;
static int num_spaces = 0;

static enum PED_Field current_field = FAMILY_ID;


#line 18 "ped_reader.c"
static const int ped_start = 1;
static const int ped_first_final = 1;
static const int ped_error = 0;

static const int ped_en_main = 1;


#line 99 "ped.ragel"



static char* get_token(char *ts, char *te)
{
    char *field = (char*) malloc ((te-ts+1) * sizeof(char));
    strncpy(field, ts, (te-ts));
    field[te-ts] = '\0';
    return field;
}

static void column_missing_error() {
    switch (current_field) {
        case FAMILY_ID:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_FAMILY_ID_MISSING));
        break;
        case INDIVIDUAL_ID:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_INDIVIDUAL_ID_MISSING));
        break;
        case FATHER_ID:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_FATHER_ID_MISSING));
        break;
        case MOTHER_ID:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_MOTHER_ID_MISSING));
        break;
        case SEX:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_SEX_MISSING));
        break;
        case PHENOTYPE:
            LOG_ERROR_F("Line %d - %s\n", lines, get_ped_syntax_error_msg(COLUMN_PHENOTYPE_MISSING));
        break;
    }

    record_error = 1;
}

static void set_field(char* ts, char *te) {
    char *field = get_token(ts, te);
    float float_val = -1.0f;
    enum Sex sex;

    if (strlen(field) == 0) {
        column_missing_error();
        return;
    }

    switch (current_field) {
        case FAMILY_ID:
            current_record = create_ped_record();
            set_ped_record_family_id(field, current_record);
        break;
        case INDIVIDUAL_ID:
            set_ped_record_individual_id(field, current_record);
            // If family is not set, assign same value as the individual ID
            if (current_record->family_id == NULL) {
                set_ped_record_family_id(field, current_record);
            }
        break;
        case FATHER_ID:
            if (strncmp("0", field, 1) != 0) {
                set_ped_record_father_id(field, current_record);
            }
        break;
        case MOTHER_ID:
            if (strncmp("0", field, 1) != 0) {
                set_ped_record_mother_id(field, current_record);
            }
        break;
        case SEX:
            sex = UNKNOWN;
            if (atoi(field) == 1) {
               sex = MALE;
            } else if (atoi(field) == 2) {
               sex = FEMALE;
            }
            set_ped_record_sex(sex, current_record);
            free(field);    // Not set as ped_record_t variable -> not freed later
        break;
        case PHENOTYPE:
            if (strncmp(".", field, 1) != 0) {
                float_val = atof(field);
            }
            set_ped_record_phenotype(float_val, current_record);
            free(field);	// Not set as ped_record_t variable -> not freed later
        break;
    }

    current_field++;
}

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

    
#line 130 "ped_reader.c"
	{
	cs = ped_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 138 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr0:
#line 96 "ped.ragel"
	{te = p+1;}
	goto st1;
tr2:
#line 72 "ped.ragel"
	{te = p+1;{
            // If batch is full, add to the list of batches and create a new, empty one
            if (ped_batch_is_full(current_batch))
            {
                list_item_t *item = list_item_new(records, 1, current_batch); 
                list_insert_item(item, batches_list);
                LOG_DEBUG_F("Batch added - %zu records", current_batch->length);
                current_batch = ped_batch_new(batch_size);
            }
            
            // If not a blank line nor erroneous record, add current record to current batch
            if (current_field != FAMILY_ID && !record_error) {
                add_record_to_ped_batch(current_record, current_batch);
                records++;
            } else {
                free(current_record);
            }
            
            lines++;
            current_field = FAMILY_ID;
            record_error = 0;
            LOG_DEBUG("End of line\n");
        }}
	goto st1;
tr7:
#line 70 "ped.ragel"
	{te = p;p--;}
	goto st1;
tr8:
#line 42 "ped.ragel"
	{te = p;p--;{
            LOG_DEBUG("Comment found, nothing to do.");
        }}
	goto st1;
tr9:
#line 46 "ped.ragel"
	{te = p;p--;{
            set_field(ts, te);
        }}
	goto st1;
tr11:
#line 1 "NONE"
	{	switch( act ) {
	case 6:
	{{p = ((te))-1;}
            set_field(ts, te);
        }
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	goto st1;
st1:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 1 "NONE"
	{ts = p;}
#line 212 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 35: goto st3;
		case 95: goto st6;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 32 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr0;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st4;
		} else
			goto st2;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st6;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st2;
			} else if ( (*p) >= 97 )
				goto st6;
		} else
			goto st2;
	} else
		goto st2;
	goto st0;
st0:
cs = 0;
	goto _out;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st2;
	goto tr7;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto tr8;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 46: goto tr10;
		case 95: goto st6;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 32 <= (*p) && (*p) <= 47 )
				goto st2;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st2;
		} else
			goto st4;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st2;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto st2;
		} else
			goto st6;
	} else
		goto st6;
	goto tr9;
tr10:
#line 1 "NONE"
	{te = p+1;}
#line 70 "ped.ragel"
	{act = 8;}
	goto st5;
tr12:
#line 1 "NONE"
	{te = p+1;}
#line 62 "ped.ragel"
	{act = 6;}
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
#line 304 "ped_reader.c"
	if ( (*p) < 48 ) {
		if ( 32 <= (*p) && (*p) <= 47 )
			goto st2;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) && (*p) <= 126 )
			goto st2;
	} else
		goto tr12;
	goto tr11;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 95 )
		goto st6;
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 32 <= (*p) && (*p) <= 47 )
				goto st2;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st2;
		} else
			goto st6;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st2;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto st2;
		} else
			goto st6;
	} else
		goto st6;
	goto tr9;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 2: goto tr7;
	case 3: goto tr8;
	case 4: goto tr9;
	case 5: goto tr11;
	case 6: goto tr9;
	}
	}

	_out: {}
	}

#line 204 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 376 "ped_reader.c"
1
#line 214 "ped.ragel"
 ) 
    {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 383 "ped_reader.c"
1
#line 217 "ped.ragel"
);
    } 

    LOG_INFO_F("Records read = %zu\n", records);

    return cs < 
#line 392 "ped_reader.c"
1
#line 222 "ped.ragel"
;
}
