
#line 1 "vcf.ragel"
#include "vcf_reader.h"

int lines = 0;
int batches = 0;
// int count = 0;


#line 11 "vcf_reader.c"
static const int vcf_start = 85;
static const int vcf_first_final = 85;
static const int vcf_error = 0;

static const int vcf_en_main = 85;


#line 305 "vcf.ragel"



int execute_vcf_ragel_machine(char *p, char *pe, list_t *batches_list, size_t batch_size, vcf_file_t *file, vcf_reader_status *status) {
    int cs;
    char *ts, *te;
    int stack[4];
    int top, act;
    char *eof = pe;

    status->current_batch->text = p;
//     printf("ragel - batch text = '%.*s'\n", 50, status->current_batch->text);

    
#line 34 "vcf_reader.c"
	{
	cs = vcf_start;
	}

#line 39 "vcf_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr4:
#line 15 "vcf.ragel"
	{
        lines++;
    }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
#line 55 "vcf_reader.c"
	switch( (*p) ) {
		case 10: goto st86;
		case 35: goto st1;
		case 95: goto tr111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr111;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr111;
	} else
		goto tr111;
	goto tr53;
tr53:
#line 104 "vcf.ragel"
	{
        printf("Line %d: Error in 'chromosome' field\n", lines);
    }
	goto st0;
tr56:
#line 120 "vcf.ragel"
	{
        printf("Line %d: Error in 'position' field\n", lines);
    }
	goto st0;
tr60:
#line 133 "vcf.ragel"
	{
        printf("Line %d: Error in 'id' field\n", lines);
    }
	goto st0;
tr64:
#line 146 "vcf.ragel"
	{
        printf("Line %d: Error in 'reference' field\n", lines);
    }
	goto st0;
tr68:
#line 163 "vcf.ragel"
	{
        printf("Line %d: Error in 'alternate' field\n", lines);
    }
	goto st0;
tr73:
#line 182 "vcf.ragel"
	{
        printf("Line %d: Error in 'quality' field\n", lines);
    }
	goto st0;
tr77:
#line 195 "vcf.ragel"
	{
        printf("Line %d: Error in 'filter' field\n", lines);
    }
	goto st0;
tr81:
#line 208 "vcf.ragel"
	{
        printf("Line %d: Error in 'info' field\n", lines);
    }
	goto st0;
tr87:
#line 221 "vcf.ragel"
	{
        printf("Line %d: Error in 'format' field\n", lines);
    }
	goto st0;
tr92:
#line 234 "vcf.ragel"
	{
        printf("Line %d: Error in sample\n", lines);
    }
	goto st0;
#line 130 "vcf_reader.c"
st0:
cs = 0;
	goto _out;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	if ( (*p) == 10 )
		goto st86;
	goto st0;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	switch( (*p) ) {
		case 35: goto st2;
		case 67: goto st4;
	}
	goto st0;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 10 )
		goto tr4;
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 72 )
		goto st5;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 82 )
		goto st6;
	goto st0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 79 )
		goto st7;
	goto st0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 77 )
		goto st8;
	goto st0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 9 )
		goto st9;
	goto st0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 80 )
		goto st10;
	goto st0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 79 )
		goto st11;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 83 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 9 )
		goto st13;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 73 )
		goto st14;
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 68 )
		goto st15;
	goto st0;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 9 )
		goto st16;
	goto st0;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( (*p) == 82 )
		goto st17;
	goto st0;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 69 )
		goto st18;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) == 70 )
		goto st19;
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 9 )
		goto st20;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 65 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 76 )
		goto st22;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 84 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 9 )
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 81 )
		goto st25;
	goto st0;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 85 )
		goto st26;
	goto st0;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 65 )
		goto st27;
	goto st0;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 76 )
		goto st28;
	goto st0;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 9 )
		goto st29;
	goto st0;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 70 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 73 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 76 )
		goto st32;
	goto st0;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( (*p) == 84 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 69 )
		goto st34;
	goto st0;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	if ( (*p) == 82 )
		goto st35;
	goto st0;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( (*p) == 9 )
		goto st36;
	goto st0;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	if ( (*p) == 73 )
		goto st37;
	goto st0;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	if ( (*p) == 78 )
		goto st38;
	goto st0;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	if ( (*p) == 70 )
		goto st39;
	goto st0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 79 )
		goto st40;
	goto st0;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	if ( (*p) == 9 )
		goto st41;
	goto st0;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	if ( (*p) == 70 )
		goto st42;
	goto st0;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	if ( (*p) == 79 )
		goto st43;
	goto st0;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	if ( (*p) == 82 )
		goto st44;
	goto st0;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	if ( (*p) == 77 )
		goto st45;
	goto st0;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	if ( (*p) == 65 )
		goto st46;
	goto st0;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 84 )
		goto st47;
	goto st0;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	if ( (*p) == 9 )
		goto st48;
	goto st0;
tr50:
#line 59 "vcf.ragel"
	{
        add_sample_name(ts, p-ts, file);
    }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 484 "vcf_reader.c"
	if ( 32 <= (*p) && (*p) <= 126 )
		goto tr49;
	goto st0;
tr49:
#line 55 "vcf.ragel"
	{
        ts = p;
    }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
#line 498 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr50;
		case 10: goto tr51;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st49;
	goto st0;
tr51:
#line 59 "vcf.ragel"
	{
        add_sample_name(ts, p-ts, file);
    }
#line 15 "vcf.ragel"
	{
        lines++;
    }
	goto st87;
tr113:
#line 230 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
#line 67 "vcf.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (vcf_batch_is_full(status->current_batch))
        {
            list_item_t *item = list_item_new(file->num_records, 1, status->current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, status->current_batch->records->size);
            status->current_batch = vcf_batch_new(batch_size);           
            if (p+1) {
                status->current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, status->current_batch->text);
            }
            batches++;
        }

        // If not a blank line, add status->current record to status->current batch
//         if (status->current_field != CHROM) {
            add_record_to_vcf_batch(status->current_record, status->current_batch);
            status->num_records++;
//             printf("Record %zu added\n", status->num_records);
//         }
        
//         status->current_field = CHROM;
        status->num_samples = 0;
//         LOG_DEBUG("\n");
    }
#line 15 "vcf.ragel"
	{
        lines++;
    }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
#line 557 "vcf_reader.c"
	switch( (*p) ) {
		case 10: goto st86;
		case 95: goto tr111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr111;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr111;
	} else
		goto tr111;
	goto tr53;
tr111:
#line 63 "vcf.ragel"
	{
        status->current_record = create_record();
    }
#line 94 "vcf.ragel"
	{
        ts = p;
//         printf("chromosome begin %c%c\n", ts[0], ts[3]);
    }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
#line 586 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr54;
		case 95: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st50;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st50;
	} else
		goto st50;
	goto tr53;
tr54:
#line 99 "vcf.ragel"
	{
//         printf("chromosome end %c%c\n", ts[0], ts[3]);
        set_record_chromosome(ts, p-ts, status->current_record);
    }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
#line 611 "vcf_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr57;
	goto tr56;
tr57:
#line 108 "vcf.ragel"
	{
        ts = p;
//         printf("position\n");
    }
	goto st52;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
#line 626 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr58;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st52;
	goto tr56;
tr58:
#line 113 "vcf.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_record_position(atol(field), status->current_record);
//         printf("Position = %ld\n", atol(field));
        free(field);
    }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
#line 645 "vcf_reader.c"
	if ( (*p) == 46 )
		goto tr61;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr62;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr62;
	} else
		goto tr62;
	goto tr60;
tr61:
#line 124 "vcf.ragel"
	{
        ts = p;
//         printf("id\n");
    }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 668 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr63;
	goto tr60;
tr63:
#line 129 "vcf.ragel"
	{
        set_record_id(ts, p-ts, status->current_record);
    }
	goto st55;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
#line 682 "vcf_reader.c"
	switch( (*p) ) {
		case 65: goto tr65;
		case 67: goto tr65;
		case 71: goto tr65;
		case 78: goto tr65;
		case 84: goto tr65;
	}
	goto tr64;
tr65:
#line 137 "vcf.ragel"
	{
        ts = p;
//         printf("reference\n");
    }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 702 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr66;
		case 65: goto st56;
		case 67: goto st56;
		case 71: goto st56;
		case 78: goto st56;
		case 84: goto st56;
	}
	goto tr64;
tr66:
#line 142 "vcf.ragel"
	{
        set_record_reference(ts, p-ts, status->current_record);
    }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 722 "vcf_reader.c"
	switch( (*p) ) {
		case 46: goto tr69;
		case 48: goto tr69;
		case 60: goto tr70;
		case 65: goto tr71;
		case 67: goto tr71;
		case 71: goto tr71;
		case 78: goto tr71;
		case 84: goto tr71;
	}
	goto tr68;
tr69:
#line 150 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
#line 745 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr72;
	goto tr68;
tr72:
#line 155 "vcf.ragel"
	{
        if (!strncmp("0", ts, 1) || !strncmp("<DEL>", ts, 5)) {
            set_record_alternate(".", 1, status->current_record);
        } else {
            set_record_alternate(ts, p-ts, status->current_record);
        }
    }
	goto st59;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
#line 763 "vcf_reader.c"
	if ( (*p) == 46 )
		goto tr74;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr75;
	goto tr73;
tr74:
#line 167 "vcf.ragel"
	{
        ts = p;
//         printf("quality\n");
    }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
#line 780 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr76;
	goto tr73;
tr76:
#line 172 "vcf.ragel"
	{
        float quality = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            quality = atof(field);
            free(field);
        }
        set_record_quality(quality, status->current_record);
    }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
#line 800 "vcf_reader.c"
	switch( (*p) ) {
		case 46: goto tr78;
		case 95: goto tr79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr79;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr79;
	} else
		goto tr79;
	goto tr77;
tr78:
#line 186 "vcf.ragel"
	{
        ts = p;
//         printf("filter\n");
    }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
#line 825 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr80;
	goto tr77;
tr80:
#line 191 "vcf.ragel"
	{
        set_record_filter(ts, p-ts, status->current_record);
    }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
#line 839 "vcf_reader.c"
	switch( (*p) ) {
		case 46: goto tr82;
		case 95: goto tr83;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr83;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr83;
	} else
		goto tr83;
	goto tr81;
tr82:
#line 199 "vcf.ragel"
	{
        ts = p;
//         printf("info\n");
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 864 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr84;
		case 59: goto st69;
		case 61: goto st71;
	}
	goto tr81;
tr84:
#line 204 "vcf.ragel"
	{
        set_record_info(ts, p-ts, status->current_record);
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 881 "vcf_reader.c"
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr88;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr88;
	} else
		goto tr88;
	goto tr87;
tr88:
#line 212 "vcf.ragel"
	{
        ts = p;
//         printf("format\n");
    }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
#line 902 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr89;
		case 58: goto st68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st66;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st66;
	} else
		goto st66;
	goto tr87;
tr89:
#line 217 "vcf.ragel"
	{
        set_record_format(ts, p-ts, status->current_record);
    }
	goto st67;
tr112:
#line 230 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 932 "vcf_reader.c"
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr93;
	goto tr92;
tr93:
#line 225 "vcf.ragel"
	{
        ts = p;
//         printf("sample\n");
    }
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
#line 947 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr112;
		case 10: goto tr113;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st88;
	goto tr92;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st66;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st66;
	} else
		goto st66;
	goto tr87;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 46: goto st64;
		case 95: goto st70;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st70;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st70;
	} else
		goto st70;
	goto tr81;
tr83:
#line 199 "vcf.ragel"
	{
        ts = p;
//         printf("info\n");
    }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
#line 996 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr84;
		case 59: goto st69;
		case 61: goto st71;
		case 95: goto st70;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st70;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st70;
	} else
		goto st70;
	goto tr81;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st72;
	goto tr81;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	if ( (*p) == 9 )
		goto tr84;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st72;
	goto tr81;
tr79:
#line 186 "vcf.ragel"
	{
        ts = p;
//         printf("filter\n");
    }
	goto st73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
#line 1039 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr80;
		case 44: goto st74;
		case 59: goto st74;
		case 95: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st73;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st73;
	} else
		goto st73;
	goto tr77;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	if ( (*p) == 95 )
		goto st73;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st73;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st73;
	} else
		goto st73;
	goto tr77;
tr75:
#line 167 "vcf.ragel"
	{
        ts = p;
//         printf("quality\n");
    }
	goto st75;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
#line 1081 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr76;
		case 46: goto st76;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st75;
	goto tr73;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st77;
	goto tr73;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	if ( (*p) == 9 )
		goto tr76;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st77;
	goto tr73;
tr70:
#line 150 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st78;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
#line 1116 "vcf_reader.c"
	if ( (*p) == 68 )
		goto st79;
	goto tr68;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	if ( (*p) == 69 )
		goto st80;
	goto tr68;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	if ( (*p) == 76 )
		goto st81;
	goto tr68;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	if ( (*p) == 62 )
		goto st58;
	goto tr68;
tr71:
#line 150 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
#line 1152 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto tr72;
		case 44: goto st83;
		case 65: goto st82;
		case 67: goto st82;
		case 71: goto st82;
		case 78: goto st82;
		case 84: goto st82;
	}
	goto tr68;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 65: goto st82;
		case 67: goto st82;
		case 71: goto st82;
		case 78: goto st82;
		case 84: goto st82;
	}
	goto tr68;
tr62:
#line 124 "vcf.ragel"
	{
        ts = p;
//         printf("id\n");
    }
	goto st84;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
#line 1186 "vcf_reader.c"
	if ( (*p) == 9 )
		goto tr63;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st84;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st84;
	} else
		goto st84;
	goto tr60;
	}
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
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
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
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
	_test_eof87: cs = 87; goto _test_eof; 
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
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 50: 
#line 104 "vcf.ragel"
	{
        printf("Line %d: Error in 'chromosome' field\n", lines);
    }
	break;
	case 51: 
	case 52: 
#line 120 "vcf.ragel"
	{
        printf("Line %d: Error in 'position' field\n", lines);
    }
	break;
	case 53: 
	case 54: 
	case 84: 
#line 133 "vcf.ragel"
	{
        printf("Line %d: Error in 'id' field\n", lines);
    }
	break;
	case 55: 
	case 56: 
#line 146 "vcf.ragel"
	{
        printf("Line %d: Error in 'reference' field\n", lines);
    }
	break;
	case 57: 
	case 58: 
	case 78: 
	case 79: 
	case 80: 
	case 81: 
	case 82: 
	case 83: 
#line 163 "vcf.ragel"
	{
        printf("Line %d: Error in 'alternate' field\n", lines);
    }
	break;
	case 59: 
	case 60: 
	case 75: 
	case 76: 
	case 77: 
#line 182 "vcf.ragel"
	{
        printf("Line %d: Error in 'quality' field\n", lines);
    }
	break;
	case 61: 
	case 62: 
	case 73: 
	case 74: 
#line 195 "vcf.ragel"
	{
        printf("Line %d: Error in 'filter' field\n", lines);
    }
	break;
	case 63: 
	case 64: 
	case 69: 
	case 70: 
	case 71: 
	case 72: 
#line 208 "vcf.ragel"
	{
        printf("Line %d: Error in 'info' field\n", lines);
    }
	break;
	case 65: 
	case 66: 
	case 68: 
#line 221 "vcf.ragel"
	{
        printf("Line %d: Error in 'format' field\n", lines);
    }
	break;
	case 67: 
#line 234 "vcf.ragel"
	{
        printf("Line %d: Error in sample\n", lines);
    }
	break;
	case 88: 
#line 230 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
#line 67 "vcf.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (vcf_batch_is_full(status->current_batch))
        {
            list_item_t *item = list_item_new(file->num_records, 1, status->current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, status->current_batch->records->size);
            status->current_batch = vcf_batch_new(batch_size);           
            if (p+1) {
                status->current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, status->current_batch->text);
            }
            batches++;
        }

        // If not a blank line, add status->current record to status->current batch
//         if (status->current_field != CHROM) {
            add_record_to_vcf_batch(status->current_record, status->current_batch);
            status->num_records++;
//             printf("Record %zu added\n", status->num_records);
//         }
        
//         status->current_field = CHROM;
        status->num_samples = 0;
//         LOG_DEBUG("\n");
    }
	break;
#line 1410 "vcf_reader.c"
	}
	}

	_out: {}
	}

#line 321 "vcf.ragel"

    
//     if (status->self_contained && !vcf_batch_is_empty(status->current_batch)) {
//     if (!vcf_batch_is_empty(status->current_batch) && !vcf_batch_is_full(status->current_batch)) {
    if (!vcf_batch_is_empty(status->current_batch)) {
        list_item_t *item = list_item_new(file->num_records, 1, status->current_batch);
        list_insert_item(item, batches_list);
//         printf("Batch added - %zu records (self-contained)\n", status->current_batch->records->size);
    }

//     printf("ragel - first chromosome = %.*s\n", 
//            ((vcf_record_t*) status->current_batch->records->items[0])->chromosome_len, 
//            ((vcf_record_t*) status->current_batch->records->items[0])->chromosome);
// 
//     printf("final state should be a minimum of %d, was %d\n",  %%{ write first_final; }%%, cs);
    return cs < 
#line 1434 "vcf_reader.c"
85
#line 336 "vcf.ragel"
;
}

vcf_reader_status *vcf_reader_status_new(size_t batch_size, int store_samples, int self_contained) {
    vcf_reader_status *status = (vcf_reader_status *) malloc (sizeof(vcf_reader_status));
    status->current_record = NULL;
    status->current_header_entry = create_header_entry();
    status->current_batch = vcf_batch_new(batch_size);

    status->num_samples = 0;
    status->num_records = 0;

    status->store_samples = store_samples;
    status->self_contained = self_contained;

    return status;
}

void vcf_reader_status_free(vcf_reader_status *status) {
    if (status->current_header_entry->name == NULL && status->current_header_entry->num_keys == 0 && status->current_header_entry->num_values == 0) {
        vcf_header_entry_free(status->current_header_entry);
    }
    if (vcf_batch_is_empty(status->current_batch)) {
        vcf_batch_free(status->current_batch);
    }
    if (status->current_record->chromosome == NULL && status->current_record->samples->size == 0) {
        vcf_record_free(status->current_record);
    }
    free(status);
}


/* **********************************************
 *              Reading and parsing             *
 * **********************************************/

int vcf_read_and_parse(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples) {
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(batch_size, read_samples, 0);
    
    if (mmap_vcf) {
        LOG_DEBUG("Using mmap for file loading\n");
        p = file->data;
        pe = p + file->data_len;
        cs = execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
    } else {
        LOG_DEBUG("Using file-IO functions for file loading\n");
        size_t max_len = 256;
        int eof_found = 0;
        char *aux;

        // Read text of a batch and call ragel parser in a loop
        while (!eof_found) {
            char *data = (char*) calloc (max_len, sizeof(char));
            int c = 0;
            int lines = 0;

            for (int i = 0; !eof_found && lines < batch_size; i++) {
                c = fgetc(file->fd);
                
                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

            }

            data[file->data_len] = '\0';

            p = data;
            pe = p + file->data_len;
            cs |= execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
            file->data_len = 0;

            // Prepare status for next batch
            status->current_batch = vcf_batch_new(batch_size);
        }
    }

    // Insert the last batch
    // TODO this should not be neccessary because it is now inserted in execute_vcf_ragel_machine
//     if (!vcf_batch_is_empty(status->current_batch))
//     {
//         list_item_t *item = list_item_new(file->num_records, 1, status->current_batch);
//         list_insert_item(item, batches_list);
//         printf("Batch added - %zu records (last)\n", status->current_batch->records->size);
//     }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", file->num_samples);

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    return cs ;
}

int vcf_gzip_read_and_parse(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples) {
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(batch_size, read_samples, 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have = 0, consumed = 0;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (consumed = 0; consumed < have && !eof_found; consumed++) {
                c = out[consumed];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    i++;
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

                // Process batch
                if (lines == batch_size) {
                    data[file->data_len] = '\0';
                    p = data;
                    pe = p + file->data_len;
                    cs |= execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
                    file->data_len = 0;

                    // Setup for next batch
                    status->current_batch = vcf_batch_new(batch_size);
                    i = 0;
                    lines = 0;
                    data = (char*) calloc (max_len, sizeof(char));
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
    if (lines > 0 && lines < batch_size) {
        data[file->data_len] = '\0';
        p = data;
        pe = p + file->data_len;
        cs |= execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
        file->data_len = 0;
    }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", file->num_samples);

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    /* clean up and return */
    (void)inflateEnd(&strm);

    return cs ;
}


/* **********************************************
 *                  Only reading                *
 * **********************************************/

int vcf_light_read(list_t *batches_list, size_t batch_size, vcf_file_t *file) {
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    char *aux;

    // Read text of a batch and call ragel parser in a loop
    while (!eof_found) {
        char *data = (char*) calloc (max_len, sizeof(char));
        int c = 0;
        int lines = 0;

        for (int i = 0; !eof_found && lines < batch_size; i++) {
            c = fgetc(file->fd);

            if (c != EOF) {
                max_len = consume_input(c, &data, max_len, i);
                if (c == '\n') {
                    lines++;
                }
            } else {
                eof_found = 1;
            }
        }

        list_item_t *item = list_item_new(file->num_records, 1, data);
        list_insert_item(item, batches_list);
//             printf("Text batch inserted = '%s'\n", data);
    }

    return 0;
}

int vcf_gzip_light_read(list_t *batches_list, size_t batch_size, vcf_file_t *file) {
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (int j = 0; j < have && !eof_found; j++) {
                c = out[j];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    i++;
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

                // Process batch
                if (lines == batch_size) {
                    list_item_t *item = list_item_new(file->num_records, 1, data);
                    list_insert_item(item, batches_list);

                    // Setup for next batch
                    i = 0;
                    lines = 0;
                    data = (char*) calloc (max_len, sizeof(char));
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
    if (lines > 0 && lines < batch_size) {
        list_item_t *item = list_item_new(file->num_records, 1, data);
        list_insert_item(item, batches_list);
    }

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}


/* **********************************************
 *      Only reading from multiple files        *
 * **********************************************/

int vcf_light_multiread(list_t **batches_list, size_t batch_size, vcf_file_t **files, size_t num_files) {
    LOG_DEBUG("Using file-IO functions for file loading\n");

    // Initialize file-private variables
    size_t max_len[num_files];
    
    for (int i = 0; i < num_files; i++) {
        max_len[i] = 256;
        files[i]->data_len = 0;
    }
    
//     char *data = NULL;
    __ssize_t line_len = 0;
    char *line = NULL;
    char *aux;
    
    int num_eof_found = 0;
    int eof_found[num_files];
    memset(eof_found, 0, num_files * sizeof(int));

    // Read text of a batch and call ragel parser in a loop
    while (num_eof_found < num_files) {
        // Read text of each file
        for (int f = 0; f < num_files; f++) {
            if (eof_found[f]) {
                printf("EOF found in file %d\n", f);
                continue;
            }

            char *data = (char*) calloc (max_len[f], sizeof(char));

            for (int i = 0; i < batch_size && !eof_found[f]; i++) {
                line_len = getline(&line, &line_len, files[f]->fd);
                if (line_len != -1) {
                    LOG_DEBUG_F("#%d Line (len %zu): %s", i, line_len, line);
                    // Line too long to be stored in data, realloc
                    if (files[f]->data_len + line_len + 1 > max_len[f]) {
                        aux = realloc(data, max_len[f] + line_len * 20);
                        if (aux) {
                            data = aux;
                            max_len[f] += line_len * 20;
                        } else {
                            LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
                        }
                    }
                    // Concat previous data with new line
                    strncat(data, line, line_len);
                    files[f]->data_len += line_len;
                } else {
                    eof_found[f] = 1;
                    num_eof_found++;
                    list_decr_writers(batches_list[f]);
                }
            }

            files[f]->data_len = 0;

            list_item_t *item = list_item_new(files[f]->num_records, 1, data);
            list_insert_item(item, batches_list[f]);
            printf("[%d] Text batch inserted\n", f);
    //             printf("Text batch inserted = '%s'\n", data);
        }
    }

    if (line != NULL) { free(line); }

    return 0;
}



/* **********************************************
 *              Auxiliary functions             *
 * **********************************************/

size_t consume_input(int c, char **data, size_t max_len, int position_in_data) {
    (*data)[position_in_data] = c;
    // Text too long to be stored in 'data', realloc
    if (position_in_data == max_len - 1) {
        char *aux = realloc(*data, max_len + 10000);
        if (aux) {
            *data = aux;
            return max_len + 10000;
        } else {
            LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
        }
    }
    return max_len;
}
