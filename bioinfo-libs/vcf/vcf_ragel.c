
#line 1 "vcf.ragel"
#include "vcf_reader.h"

int lines = 0;
int batches = 0;
// int count = 0;


#line 11 "vcf_ragel.c"
static const int vcf_start = 100;
static const int vcf_first_final = 100;
static const int vcf_error = 0;

static const int vcf_en_main = 100;


#line 301 "vcf.ragel"



int execute_vcf_ragel_machine(char *p, char *pe, list_t *batches_list, size_t batch_size, vcf_file_t *file, vcf_reader_status *status) {
    int cs;
    char *ts, *te;
    int stack[4];
    int top, act;
    char *eof = pe;

    status->current_batch->text = p;
//     printf("ragel - batch text = '%.*s'\n", 50, status->current_batch->text);

    
#line 34 "vcf_ragel.c"
	{
	cs = vcf_start;
	}

#line 39 "vcf_ragel.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 100:
	switch( (*p) ) {
		case 10: goto st101;
		case 35: goto st86;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
tr53:
#line 105 "vcf.ragel"
	{
        printf("Line %d: Error in 'chromosome' field\n", lines);
    }
	goto st0;
tr56:
#line 121 "vcf.ragel"
	{
        printf("Line %d: Error in 'position' field\n", lines);
    }
	goto st0;
tr60:
#line 134 "vcf.ragel"
	{
        printf("Line %d: Error in 'id' field\n", lines);
    }
	goto st0;
tr64:
#line 147 "vcf.ragel"
	{
        printf("Line %d: Error in 'reference' field\n", lines);
    }
	goto st0;
tr68:
#line 164 "vcf.ragel"
	{
        printf("Line %d: Error in 'alternate' field\n", lines);
    }
	goto st0;
tr73:
#line 183 "vcf.ragel"
	{
        printf("Line %d: Error in 'quality' field\n", lines);
    }
	goto st0;
tr77:
#line 196 "vcf.ragel"
	{
        printf("Line %d: Error in 'filter' field\n", lines);
    }
	goto st0;
tr81:
#line 209 "vcf.ragel"
	{
        printf("Line %d: Error in 'info' field\n", lines);
    }
	goto st0;
tr87:
#line 222 "vcf.ragel"
	{
        printf("Line %d: Error in 'format' field\n", lines);
    }
	goto st0;
tr92:
#line 235 "vcf.ragel"
	{
        printf("Line %d: Error in sample\n", lines);
    }
	goto st0;
#line 120 "vcf_ragel.c"
st0:
cs = 0;
	goto _out;
tr122:
#line 15 "vcf.ragel"
	{
        lines++;
    }
	goto st101;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
#line 134 "vcf_ragel.c"
	switch( (*p) ) {
		case 10: goto st101;
		case 35: goto st1;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	switch( (*p) ) {
		case 35: goto st2;
		case 67: goto st5;
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
tr4:
#line 15 "vcf.ragel"
	{
        lines++;
    }
	goto st102;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
#line 184 "vcf_ragel.c"
	switch( (*p) ) {
		case 10: goto st103;
		case 35: goto st1;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 10: goto st103;
		case 35: goto st4;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 67 )
		goto st5;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 72 )
		goto st6;
	goto st0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 82 )
		goto st7;
	goto st0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 79 )
		goto st8;
	goto st0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 77 )
		goto st9;
	goto st0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 9 )
		goto st10;
	goto st0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 80 )
		goto st11;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 79 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 83 )
		goto st13;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 9 )
		goto st14;
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 73 )
		goto st15;
	goto st0;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 68 )
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
	if ( (*p) == 82 )
		goto st18;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) == 69 )
		goto st19;
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 70 )
		goto st20;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 9 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 65 )
		goto st22;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 76 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 84 )
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 9 )
		goto st25;
	goto st0;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 81 )
		goto st26;
	goto st0;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 85 )
		goto st27;
	goto st0;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 65 )
		goto st28;
	goto st0;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 76 )
		goto st29;
	goto st0;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 9 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 70 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 73 )
		goto st32;
	goto st0;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( (*p) == 76 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 84 )
		goto st34;
	goto st0;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	if ( (*p) == 69 )
		goto st35;
	goto st0;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( (*p) == 82 )
		goto st36;
	goto st0;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	if ( (*p) == 9 )
		goto st37;
	goto st0;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	if ( (*p) == 73 )
		goto st38;
	goto st0;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	if ( (*p) == 78 )
		goto st39;
	goto st0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 70 )
		goto st40;
	goto st0;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	if ( (*p) == 79 )
		goto st41;
	goto st0;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	if ( (*p) == 9 )
		goto st42;
	goto st0;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	if ( (*p) == 70 )
		goto st43;
	goto st0;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	if ( (*p) == 79 )
		goto st44;
	goto st0;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	if ( (*p) == 82 )
		goto st45;
	goto st0;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	if ( (*p) == 77 )
		goto st46;
	goto st0;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 65 )
		goto st47;
	goto st0;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	if ( (*p) == 84 )
		goto st48;
	goto st0;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	if ( (*p) == 9 )
		goto st49;
	goto st0;
tr50:
#line 59 "vcf.ragel"
	{
        add_sample_name(ts, p-ts, file);
    }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
#line 542 "vcf_ragel.c"
	if ( 32 <= (*p) && (*p) <= 126 )
		goto tr49;
	goto st0;
tr49:
#line 55 "vcf.ragel"
	{
        ts = p;
    }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
#line 556 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr50;
		case 10: goto tr51;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st50;
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
	goto st104;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
#line 578 "vcf_ragel.c"
	switch( (*p) ) {
		case 10: goto st104;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
tr125:
#line 63 "vcf.ragel"
	{
        status->current_record = create_record();
    }
#line 95 "vcf.ragel"
	{
        ts = p;
//         printf("chromosome begin %c%c\n", ts[0], ts[3]);
    }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
#line 607 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr54;
		case 95: goto st51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st51;
	} else
		goto st51;
	goto tr53;
tr54:
#line 100 "vcf.ragel"
	{
//         printf("chromosome end %c...%c\n", *ts, *p);
        set_record_chromosome(ts, p-ts, status->current_record);
    }
	goto st52;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
#line 632 "vcf_ragel.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr57;
	goto tr56;
tr57:
#line 109 "vcf.ragel"
	{
        ts = p;
//         printf("position\n");
    }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
#line 647 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr58;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st53;
	goto tr56;
tr58:
#line 114 "vcf.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_record_position(atol(field), status->current_record);
//         printf("Position = %ld\n", atol(field));
        free(field);
    }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 666 "vcf_ragel.c"
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
#line 125 "vcf.ragel"
	{
        ts = p;
//         printf("id\n");
    }
	goto st55;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
#line 689 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr63;
	goto tr60;
tr63:
#line 130 "vcf.ragel"
	{
        set_record_id(ts, p-ts, status->current_record);
    }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 703 "vcf_ragel.c"
	switch( (*p) ) {
		case 65: goto tr65;
		case 67: goto tr65;
		case 71: goto tr65;
		case 78: goto tr65;
		case 84: goto tr65;
	}
	goto tr64;
tr65:
#line 138 "vcf.ragel"
	{
        ts = p;
//         printf("reference\n");
    }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 723 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr66;
		case 65: goto st57;
		case 67: goto st57;
		case 71: goto st57;
		case 78: goto st57;
		case 84: goto st57;
	}
	goto tr64;
tr66:
#line 143 "vcf.ragel"
	{
        set_record_reference(ts, p-ts, status->current_record);
    }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
#line 743 "vcf_ragel.c"
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
#line 151 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st59;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
#line 766 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr72;
	goto tr68;
tr72:
#line 156 "vcf.ragel"
	{
        if (!strncmp("0", ts, 1) || !strncmp("<DEL>", ts, 5)) {
            set_record_alternate(".", 1, status->current_record);
        } else {
            set_record_alternate(ts, p-ts, status->current_record);
        }
    }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
#line 784 "vcf_ragel.c"
	if ( (*p) == 46 )
		goto tr74;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr75;
	goto tr73;
tr74:
#line 168 "vcf.ragel"
	{
        ts = p;
//         printf("quality\n");
    }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
#line 801 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr76;
	goto tr73;
tr76:
#line 173 "vcf.ragel"
	{
        float quality = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            quality = atof(field);
            free(field);
        }
        set_record_quality(quality, status->current_record);
    }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
#line 821 "vcf_ragel.c"
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
#line 187 "vcf.ragel"
	{
        ts = p;
//         printf("filter\n");
    }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
#line 846 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr80;
	goto tr77;
tr80:
#line 192 "vcf.ragel"
	{
        set_record_filter(ts, p-ts, status->current_record);
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 860 "vcf_ragel.c"
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
#line 200 "vcf.ragel"
	{
        ts = p;
//         printf("info\n");
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 885 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr84;
		case 59: goto st70;
		case 61: goto st72;
	}
	goto tr81;
tr84:
#line 205 "vcf.ragel"
	{
        set_record_info(ts, p-ts, status->current_record);
    }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
#line 902 "vcf_ragel.c"
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
#line 213 "vcf.ragel"
	{
        ts = p;
//         printf("format\n");
    }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 923 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr89;
		case 58: goto st69;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st67;
	} else
		goto st67;
	goto tr87;
tr89:
#line 218 "vcf.ragel"
	{
        set_record_format(ts, p-ts, status->current_record);
    }
	goto st68;
tr130:
#line 231 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
#line 953 "vcf_ragel.c"
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr93;
	goto tr92;
tr93:
#line 226 "vcf.ragel"
	{
        ts = p;
//         printf("sample\n");
    }
	goto st105;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
#line 968 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr130;
		case 10: goto tr131;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st105;
	goto tr92;
tr131:
#line 231 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
#line 67 "vcf.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
//         if (vcf_batch_is_full(status->current_batch))
        if (batch_size > 0 && status->current_batch->records->size == batch_size)
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
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
#line 1018 "vcf_ragel.c"
	switch( (*p) ) {
		case 10: goto st107;
		case 95: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr125;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr125;
	} else
		goto tr125;
	goto tr53;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	if ( (*p) == 10 )
		goto st107;
	goto st0;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st67;
	} else
		goto st67;
	goto tr87;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 46: goto st65;
		case 95: goto st71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st71;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st71;
	} else
		goto st71;
	goto tr81;
tr83:
#line 200 "vcf.ragel"
	{
        ts = p;
//         printf("info\n");
    }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
#line 1080 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr84;
		case 59: goto st70;
		case 61: goto st72;
		case 95: goto st71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st71;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st71;
	} else
		goto st71;
	goto tr81;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st73;
	goto tr81;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	if ( (*p) == 9 )
		goto tr84;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st73;
	goto tr81;
tr79:
#line 187 "vcf.ragel"
	{
        ts = p;
//         printf("filter\n");
    }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
#line 1123 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr80;
		case 44: goto st75;
		case 59: goto st75;
		case 95: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st74;
	} else
		goto st74;
	goto tr77;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	if ( (*p) == 95 )
		goto st74;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st74;
	} else
		goto st74;
	goto tr77;
tr75:
#line 168 "vcf.ragel"
	{
        ts = p;
//         printf("quality\n");
    }
	goto st76;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
#line 1165 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr76;
		case 46: goto st77;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st76;
	goto tr73;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st78;
	goto tr73;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	if ( (*p) == 9 )
		goto tr76;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st78;
	goto tr73;
tr70:
#line 151 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
#line 1200 "vcf_ragel.c"
	if ( (*p) == 68 )
		goto st80;
	goto tr68;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	if ( (*p) == 69 )
		goto st81;
	goto tr68;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	if ( (*p) == 76 )
		goto st82;
	goto tr68;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	if ( (*p) == 62 )
		goto st59;
	goto tr68;
tr71:
#line 151 "vcf.ragel"
	{
        ts = p;
//         printf("alternate\n");
    }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
#line 1236 "vcf_ragel.c"
	switch( (*p) ) {
		case 9: goto tr72;
		case 44: goto st84;
		case 65: goto st83;
		case 67: goto st83;
		case 71: goto st83;
		case 78: goto st83;
		case 84: goto st83;
	}
	goto tr68;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 65: goto st83;
		case 67: goto st83;
		case 71: goto st83;
		case 78: goto st83;
		case 84: goto st83;
	}
	goto tr68;
tr62:
#line 125 "vcf.ragel"
	{
        ts = p;
//         printf("id\n");
    }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
#line 1270 "vcf_ragel.c"
	if ( (*p) == 9 )
		goto tr63;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st85;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st85;
	} else
		goto st85;
	goto tr60;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 35: goto st87;
		case 67: goto st5;
	}
	goto st0;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	if ( (*p) == 102 )
		goto st88;
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 10: goto tr4;
		case 105: goto st89;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 10: goto tr4;
		case 108: goto st90;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 10: goto tr4;
		case 101: goto st91;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 10: goto tr4;
		case 102: goto st92;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 10: goto tr4;
		case 111: goto st93;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 10: goto tr4;
		case 114: goto st94;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 10: goto tr4;
		case 109: goto st95;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	switch( (*p) ) {
		case 10: goto tr4;
		case 97: goto st96;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 10: goto tr4;
		case 116: goto st97;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 10: goto tr4;
		case 61: goto st98;
	}
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st3;
	goto st0;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 10: goto tr4;
		case 32: goto st3;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st99;
	goto st0;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 10: goto tr122;
		case 32: goto st3;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st99;
	goto st0;
	}
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
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
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
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
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
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
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 51: 
#line 105 "vcf.ragel"
	{
        printf("Line %d: Error in 'chromosome' field\n", lines);
    }
	break;
	case 52: 
	case 53: 
#line 121 "vcf.ragel"
	{
        printf("Line %d: Error in 'position' field\n", lines);
    }
	break;
	case 54: 
	case 55: 
	case 85: 
#line 134 "vcf.ragel"
	{
        printf("Line %d: Error in 'id' field\n", lines);
    }
	break;
	case 56: 
	case 57: 
#line 147 "vcf.ragel"
	{
        printf("Line %d: Error in 'reference' field\n", lines);
    }
	break;
	case 58: 
	case 59: 
	case 79: 
	case 80: 
	case 81: 
	case 82: 
	case 83: 
	case 84: 
#line 164 "vcf.ragel"
	{
        printf("Line %d: Error in 'alternate' field\n", lines);
    }
	break;
	case 60: 
	case 61: 
	case 76: 
	case 77: 
	case 78: 
#line 183 "vcf.ragel"
	{
        printf("Line %d: Error in 'quality' field\n", lines);
    }
	break;
	case 62: 
	case 63: 
	case 74: 
	case 75: 
#line 196 "vcf.ragel"
	{
        printf("Line %d: Error in 'filter' field\n", lines);
    }
	break;
	case 64: 
	case 65: 
	case 70: 
	case 71: 
	case 72: 
	case 73: 
#line 209 "vcf.ragel"
	{
        printf("Line %d: Error in 'info' field\n", lines);
    }
	break;
	case 66: 
	case 67: 
	case 69: 
#line 222 "vcf.ragel"
	{
        printf("Line %d: Error in 'format' field\n", lines);
    }
	break;
	case 68: 
#line 235 "vcf.ragel"
	{
        printf("Line %d: Error in sample\n", lines);
    }
	break;
	case 105: 
#line 231 "vcf.ragel"
	{
        add_record_sample(ts, p-ts, status->current_record);
    }
#line 67 "vcf.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
//         if (vcf_batch_is_full(status->current_batch))
        if (batch_size > 0 && status->current_batch->records->size == batch_size)
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
#line 1663 "vcf_ragel.c"
	}
	}

	_out: {}
	}

#line 317 "vcf.ragel"

    
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
#line 1685 "vcf_ragel.c"
100
#line 330 "vcf.ragel"
;
}
