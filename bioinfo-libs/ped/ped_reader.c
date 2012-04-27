
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


#line 107 "ped.ragel"



static char* get_token(char *ts, char *te)
{
    char *field = (char*) malloc ((te-ts+1) * sizeof(char));
    strncpy(field, ts, (te-ts));
    field[te-ts] = '\0';
    return field;
}

static void set_field(char* ts, char *te) {
    char *field = get_token(ts, te);
    float float_val = -1.0f;
    enum Sex sex;

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

//     LOG_DEBUG_F("Line %zu, current_field = %d\n", lines, current_field);
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

    
#line 101 "ped_reader.c"
	{
	cs = ped_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 109 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr0:
#line 104 "ped.ragel"
	{te = p+1;}
	goto st1;
tr2:
#line 77 "ped.ragel"
	{te = p+1;{
            // If batch is full, add to the list of batches and create a new, empty one
            if (ped_batch_is_full(current_batch))
            {
                list_item_t *item = list_item_new(records, 1, current_batch);
                list_insert_item(item, batches_list);
                LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
                current_batch = ped_batch_new(batch_size);
            }
            
            // If not a blank line nor erroneous record, add current record to current batch
            if (current_field > PHENOTYPE) {
                add_record_to_ped_batch(current_record, current_batch);
                records++;
            } else {
                if (current_field > FAMILY_ID) {
                    LOG_ERROR_F("Line %d - %d field(s) missing\n", lines, OTHER - current_field);
                    free(current_record);
                }
            }
            
            lines++;
            current_field = FAMILY_ID;
            record_error = 0;
            LOG_DEBUG("End of line\n");
        }}
	goto st1;
tr8:
#line 49 "ped.ragel"
	{te = p;p--;{
            set_field(ts, te);
        }}
	goto st1;
tr9:
#line 45 "ped.ragel"
	{te = p;p--;{
            LOG_DEBUG("Comment found, nothing to do.");
        }}
	goto st1;
tr11:
#line 73 "ped.ragel"
	{te = p;p--;{
            LOG_DEBUG("Genotypes found, nothing to do.");
        }}
	goto st1;
st1:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 1 "NONE"
	{ts = p;}
#line 174 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 32: goto tr0;
		case 35: goto st3;
		case 66: goto st6;
		case 71: goto st7;
		case 78: goto st7;
		case 84: goto st7;
		case 95: goto st6;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr0;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto st2;
			} else if ( (*p) >= 48 )
				goto st4;
		} else
			goto st2;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
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
		goto st7;
	goto st0;
st0:
cs = 0;
	goto _out;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st2;
	goto tr8;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto tr9;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 46: goto st5;
		case 95: goto st6;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
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
	goto tr8;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) < 48 ) {
		if ( 33 <= (*p) && (*p) <= 47 )
			goto st2;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) && (*p) <= 126 )
			goto st2;
	} else
		goto st5;
	goto tr8;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 95 )
		goto st6;
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
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
	goto tr8;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	switch( (*p) ) {
		case 65: goto st7;
		case 67: goto st7;
		case 71: goto st7;
		case 78: goto st7;
		case 84: goto st7;
	}
	goto tr11;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 2: goto tr8;
	case 3: goto tr9;
	case 4: goto tr8;
	case 5: goto tr8;
	case 6: goto tr8;
	case 7: goto tr11;
	}
	}

	_out: {}
	}

#line 183 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 347 "ped_reader.c"
1
#line 193 "ped.ragel"
 ) 
    {
        LOG_ERROR("The file was not successfully read\n");
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 355 "ped_reader.c"
1
#line 197 "ped.ragel"
);
    } 

    LOG_INFO_F("Records read = %zu\n", records);

    return cs < 
#line 364 "ped_reader.c"
1
#line 202 "ped.ragel"
;
}
