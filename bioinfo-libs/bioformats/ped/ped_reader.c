
#line 1 "ped.ragel"
#include "ped_reader.h"

static size_t lines = 0;
static size_t num_records = 0;
static size_t genotype = 0;

static ped_record_t *current_record;
static ped_batch_t *current_batch;


#line 14 "ped_reader.c"
static const int ped_start = 1;
static const int ped_first_final = 21;
static const int ped_error = 0;

static const int ped_en_main = 1;


#line 205 "ped.ragel"




int ped_ragel_read(list_t *batches_list, size_t batch_size, ped_file_t *file)
{
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;
    int custom_field_count = 0;

    current_batch = ped_batch_new(batch_size);

    
#line 41 "ped_reader.c"
	{
	cs = ped_start;
	}

#line 46 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 1:
	if ( (*p) == 35 )
		goto st17;
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 10 )
		goto st22;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr37;
	goto tr2;
tr2:
#line 53 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	goto st0;
tr5:
#line 65 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	goto st0;
tr9:
#line 79 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr13:
#line 93 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr17:
#line 113 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr21:
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
tr28:
#line 145 "ped.ragel"
	{
        printf("Line %zu: Error in 'header' field\n", lines);
    }
	goto st0;
tr42:
#line 165 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
#line 113 "ped_reader.c"
st0:
cs = 0;
	goto _out;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 10 )
		goto st22;
	goto st0;
tr37:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 45 "ped.ragel"
	{
        ts = p;
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 139 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr3;
		case 32: goto tr3;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st2;
	goto tr2;
tr3:
#line 49 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
#line 157 "ped_reader.c"
	if ( (*p) == 95 )
		goto tr6;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr6;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr6;
	} else
		goto tr6;
	goto tr5;
tr6:
#line 57 "ped.ragel"
	{
        ts = p;
    }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 179 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr7;
		case 32: goto tr7;
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
	goto tr5;
tr7:
#line 61 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
#line 204 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr10;
		case 95: goto tr11;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr11;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr11;
	} else
		goto tr11;
	goto tr9;
tr10:
#line 69 "ped.ragel"
	{
        ts = p;
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 228 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr12;
		case 32: goto tr12;
	}
	goto tr9;
tr12:
#line 73 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 246 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr14;
		case 95: goto tr15;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr15;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr15;
	} else
		goto tr15;
	goto tr13;
tr14:
#line 83 "ped.ragel"
	{
        ts = p;
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 270 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr16;
		case 32: goto tr16;
	}
	goto tr13;
tr16:
#line 87 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 288 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr18;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr19;
	goto tr17;
tr18:
#line 97 "ped.ragel"
	{
        ts = p;
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 304 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr20;
		case 32: goto tr20;
	}
	goto tr17;
tr20:
#line 101 "ped.ragel"
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
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 328 "ped_reader.c"
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr22;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr22;
	} else
		goto tr22;
	goto tr21;
tr22:
#line 117 "ped.ragel"
	{
        ts = p;
    }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 348 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr38;
		case 10: goto tr39;
		case 32: goto tr38;
	}
	if ( (*p) < 48 ) {
		if ( 11 <= (*p) && (*p) <= 13 )
			goto tr40;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st23;
		} else if ( (*p) >= 65 )
			goto st23;
	} else
		goto st23;
	goto tr21;
tr38:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 149 "ped.ragel"
	{
        custom_field_count = 6;
    }
	goto st24;
tr46:
#line 157 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
#line 393 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr44;
		case 32: goto st25;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto tr45;
	} else if ( (*p) >= 9 )
		goto st25;
	goto tr42;
tr40:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 149 "ped.ragel"
	{
        custom_field_count = 6;
    }
	goto st25;
tr48:
#line 157 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 431 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr44;
		case 32: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st25;
	goto st0;
tr39:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 149 "ped.ragel"
	{
        custom_field_count = 6;
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st26;
tr44:
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st26;
tr47:
#line 157 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 534 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr44;
		case 32: goto st25;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto tr37;
	} else if ( (*p) >= 9 )
		goto st25;
	goto tr2;
tr45:
#line 153 "ped.ragel"
	{
        ts = p;
    }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 555 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr46;
		case 10: goto tr47;
		case 32: goto tr46;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto st27;
	} else if ( (*p) >= 11 )
		goto tr48;
	goto tr42;
tr19:
#line 97 "ped.ragel"
	{
        ts = p;
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
#line 577 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr20;
		case 32: goto tr20;
		case 46: goto st13;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st12;
	goto tr17;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st14;
	goto tr17;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	switch( (*p) ) {
		case 9: goto tr20;
		case 32: goto tr20;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st14;
	goto tr17;
tr15:
#line 83 "ped.ragel"
	{
        ts = p;
    }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 614 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr16;
		case 32: goto tr16;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st15;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st15;
	} else
		goto st15;
	goto tr13;
tr11:
#line 69 "ped.ragel"
	{
        ts = p;
    }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 639 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr12;
		case 32: goto tr12;
		case 95: goto st16;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st16;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st16;
	} else
		goto st16;
	goto tr9;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 9: goto st18;
		case 10: goto st19;
		case 32: goto st18;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr31;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr31;
	} else
		goto tr31;
	goto tr28;
tr33:
#line 136 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (!strncmp(field_name, file->custom_field,sizeof(field_name))) {
            file->num_field = custom_field_count;
        }
        free(field_name);
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 687 "ped_reader.c"
	if ( (*p) == 10 )
		goto st19;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr31;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr31;
	} else
		goto tr31;
	goto tr28;
tr34:
#line 136 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (!strncmp(field_name, file->custom_field,sizeof(field_name))) {
            file->num_field = custom_field_count;
        }
        free(field_name);
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 714 "ped_reader.c"
	if ( (*p) == 35 )
		goto st0;
	goto st21;
tr31:
#line 132 "ped.ragel"
	{
        ts = p;
    }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 728 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr33;
		case 10: goto tr34;
		case 32: goto tr33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st20;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st20;
	} else
		goto st20;
	goto tr28;
	}
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
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
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 24: 
	case 25: 
	case 26: 
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
	case 2: 
#line 53 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	break;
	case 3: 
	case 4: 
#line 65 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	break;
	case 5: 
	case 6: 
	case 16: 
#line 79 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 7: 
	case 8: 
	case 15: 
#line 93 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 9: 
	case 10: 
	case 12: 
	case 13: 
	case 14: 
#line 113 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 11: 
#line 128 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 17: 
	case 18: 
	case 20: 
#line 145 "ped.ragel"
	{
        printf("Line %zu: Error in 'header' field\n", lines);
    }
	break;
	case 27: 
#line 157 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
	case 23: 
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 149 "ped.ragel"
	{
        custom_field_count = 6;
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
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
#line 909 "ped_reader.c"
	}
	}

	_out: {}
	}

#line 225 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 928 "ped_reader.c"
21
#line 235 "ped.ragel"
 ) 
    {
        LOG_ERROR("The file was not successfully read\n");
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 936 "ped_reader.c"
21
#line 239 "ped.ragel"
);
    } 

    LOG_INFO_F("PED records read = %zu\n", num_records);

    return cs < 
#line 945 "ped_reader.c"
21
#line 244 "ped.ragel"
;
}
