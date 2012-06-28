
#line 1 "vcf.ragel"
#include "vcf_reader.h"


#line 7 "vcf_reader.c"
static const int vcf_start = 82;
static const int vcf_first_final = 82;
static const int vcf_error = 0;

static const int vcf_en_fileformat_scan = 86;
static const int vcf_en_header_scan = 97;
static const int vcf_en_delimiter_scan = 104;
static const int vcf_en_record_scan = 106;
static const int vcf_en_main = 82;


#line 230 "vcf.ragel"


static char* get_token(char *ts, char *te) {
    return strndup(ts, te-ts);
}

static void set_field(char* ts, char *te, vcf_reader_status *status) {
    char *field = get_token(ts, te);
    float quality = -1.0f;

    switch (status->current_field) {
        case CHROM:
            status->current_record = create_record();
            set_record_chromosome(field, status->current_record);
        break;
        case POS:
            set_record_position(atol(field), status->current_record);
            free(field);	// Not set as vcf_record_t variable -> not freed later
        break;
        case ID:
            set_record_id(field, status->current_record);
        break;
        case REF:
            set_record_reference(field, status->current_record);
        break;
        case ALT:
            if (!strcmp("0", field) || !strcmp("<DEL>", field)) {
                free(field);
                field = strdup(".");
            }
            set_record_alternate(field, status->current_record);
        break;
        case QUAL:
            if (strncmp(".", field, 1) != 0) {
                quality = atof(field);
            }
            set_record_quality(quality, status->current_record);
            free(field);	// Not set as vcf_record_t variable -> not freed later
        break;
        case FILTER:
            set_record_filter(field, status->current_record);
        break;
        case INFO:
            set_record_info(field, status->current_record);
        break;
        case FORMAT:
            set_record_format(field, status->current_record);
        break;
        case SAMPLE:
            if(status->store_samples) { 
                add_record_sample(field, status->current_record, &(status->num_samples)); 
            } else {
                free(field);    // Not added to samples -> not freed later
            }
        break;
    }

    if (status->current_field < SAMPLE) {
        status->current_field++;
    }
}

int execute_vcf_ragel_machine(char *p, char *pe, list_t *batches_list, size_t batch_size, vcf_file_t *file, vcf_reader_status *status) {
    int cs;
    char *ts, *te;
    int stack[4];
    int top, act;
    char *eof = pe;

    
#line 90 "vcf_reader.c"
	{
	cs = vcf_start;
	top = 0;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 99 "vcf_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	goto _resume;

_again:
	switch ( cs ) {
		case 82: goto st82;
		case 0: goto st0;
		case 1: goto st1;
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
		case 18: goto st18;
		case 19: goto st19;
		case 20: goto st20;
		case 21: goto st21;
		case 22: goto st22;
		case 23: goto st23;
		case 24: goto st24;
		case 25: goto st25;
		case 26: goto st26;
		case 27: goto st27;
		case 28: goto st28;
		case 29: goto st29;
		case 30: goto st30;
		case 31: goto st31;
		case 32: goto st32;
		case 33: goto st33;
		case 34: goto st34;
		case 35: goto st35;
		case 36: goto st36;
		case 37: goto st37;
		case 38: goto st38;
		case 39: goto st39;
		case 40: goto st40;
		case 41: goto st41;
		case 42: goto st42;
		case 43: goto st43;
		case 44: goto st44;
		case 45: goto st45;
		case 46: goto st46;
		case 47: goto st47;
		case 48: goto st48;
		case 83: goto st83;
		case 49: goto st49;
		case 50: goto st50;
		case 51: goto st51;
		case 52: goto st52;
		case 53: goto st53;
		case 54: goto st54;
		case 55: goto st55;
		case 56: goto st56;
		case 57: goto st57;
		case 58: goto st58;
		case 59: goto st59;
		case 60: goto st60;
		case 61: goto st61;
		case 62: goto st62;
		case 63: goto st63;
		case 64: goto st64;
		case 65: goto st65;
		case 66: goto st66;
		case 84: goto st84;
		case 85: goto st85;
		case 67: goto st67;
		case 68: goto st68;
		case 69: goto st69;
		case 70: goto st70;
		case 71: goto st71;
		case 72: goto st72;
		case 73: goto st73;
		case 74: goto st74;
		case 75: goto st75;
		case 76: goto st76;
		case 77: goto st77;
		case 78: goto st78;
		case 79: goto st79;
		case 80: goto st80;
		case 86: goto st86;
		case 87: goto st87;
		case 88: goto st88;
		case 89: goto st89;
		case 90: goto st90;
		case 91: goto st91;
		case 92: goto st92;
		case 93: goto st93;
		case 94: goto st94;
		case 95: goto st95;
		case 96: goto st96;
		case 97: goto st97;
		case 98: goto st98;
		case 99: goto st99;
		case 81: goto st81;
		case 100: goto st100;
		case 101: goto st101;
		case 102: goto st102;
		case 103: goto st103;
		case 104: goto st104;
		case 105: goto st105;
		case 106: goto st106;
		case 107: goto st107;
		case 108: goto st108;
		case 109: goto st109;
		case 110: goto st110;
		case 111: goto st111;
		case 112: goto st112;
		case 113: goto st113;
		case 114: goto st114;
		case 115: goto st115;
		case 116: goto st116;
		case 117: goto st117;
		case 118: goto st118;
		case 119: goto st119;
		case 120: goto st120;
		case 121: goto st121;
		case 122: goto st122;
		case 123: goto st123;
		case 124: goto st124;
		case 125: goto st125;
	default: break;
	}

	if ( ++p == pe )
		goto _test_eof;
_resume:
	switch ( cs )
	{
st82:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof82;
case 82:
#line 247 "vcf_reader.c"
	switch( (*p) ) {
		case 35: goto st1;
		case 95: goto tr86;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr86;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr86;
	} else
		goto tr86;
	goto st0;
st0:
cs = 0;
	goto _out;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 35 )
		goto tr0;
	goto st0;
tr0:
#line 20 "vcf.ragel"
	{ {stack[top++] = 2; goto st86;} }
	goto st2;
st2:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 281 "vcf_reader.c"
	if ( (*p) == 35 )
		goto st3;
	goto st0;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 35 )
		goto tr3;
	goto st0;
tr3:
#line 50 "vcf.ragel"
	{ {stack[top++] = 4; goto st97;} }
	goto st4;
st4:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 302 "vcf_reader.c"
	if ( (*p) == 35 )
		goto st5;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 35: goto tr3;
		case 67: goto st6;
	}
	goto st0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 72 )
		goto st7;
	goto st0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 82 )
		goto st8;
	goto st0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 79 )
		goto st9;
	goto st0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 77 )
		goto st10;
	goto st0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 9 )
		goto st11;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 80 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 79 )
		goto st13;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 83 )
		goto st14;
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
	if ( (*p) == 73 )
		goto st16;
	goto st0;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( (*p) == 68 )
		goto st17;
	goto st0;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 9 )
		goto st18;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) == 82 )
		goto st19;
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 69 )
		goto st20;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 70 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 9 )
		goto st22;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 65 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 76 )
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 84 )
		goto st25;
	goto st0;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 9 )
		goto st26;
	goto st0;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 81 )
		goto st27;
	goto st0;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 85 )
		goto st28;
	goto st0;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 65 )
		goto st29;
	goto st0;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 76 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 9 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 70 )
		goto st32;
	goto st0;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( (*p) == 73 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 76 )
		goto st34;
	goto st0;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	if ( (*p) == 84 )
		goto st35;
	goto st0;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( (*p) == 69 )
		goto st36;
	goto st0;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	if ( (*p) == 82 )
		goto st37;
	goto st0;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	if ( (*p) == 9 )
		goto st38;
	goto st0;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	if ( (*p) == 73 )
		goto st39;
	goto st0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 78 )
		goto st40;
	goto st0;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	if ( (*p) == 70 )
		goto st41;
	goto st0;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	if ( (*p) == 79 )
		goto st42;
	goto st0;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	if ( (*p) == 9 )
		goto st43;
	goto st0;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	if ( (*p) == 70 )
		goto st44;
	goto st0;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	if ( (*p) == 79 )
		goto st45;
	goto st0;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	if ( (*p) == 82 )
		goto st46;
	goto st0;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 77 )
		goto st47;
	goto st0;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	if ( (*p) == 65 )
		goto st48;
	goto st0;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	if ( (*p) == 84 )
		goto tr48;
	goto st0;
tr48:
#line 111 "vcf.ragel"
	{ {stack[top++] = 83; goto st104;} }
	goto st83;
st83:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof83;
case 83:
#line 626 "vcf_reader.c"
	if ( (*p) == 95 )
		goto tr86;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr86;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr86;
	} else
		goto tr86;
	goto st0;
tr86:
#line 227 "vcf.ragel"
	{ p--; {stack[top++] = 49; goto st106;} }
	goto st49;
st49:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof49;
case 49:
#line 648 "vcf_reader.c"
	switch( (*p) ) {
		case 9: goto st50;
		case 95: goto st49;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st49;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st49;
	} else
		goto st49;
	goto st0;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st51;
	goto st0;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	if ( (*p) == 9 )
		goto st52;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st51;
	goto st0;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	if ( (*p) == 46 )
		goto st53;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st80;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st80;
	} else
		goto st80;
	goto st0;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	if ( (*p) == 9 )
		goto st54;
	goto st0;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 65: goto st55;
		case 67: goto st55;
		case 71: goto st55;
		case 78: goto st55;
		case 84: goto st55;
	}
	goto st0;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 9: goto st56;
		case 65: goto st55;
		case 67: goto st55;
		case 71: goto st55;
		case 78: goto st55;
		case 84: goto st55;
	}
	goto st0;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 46: goto st57;
		case 48: goto st57;
		case 60: goto st74;
		case 65: goto st78;
		case 67: goto st78;
		case 71: goto st78;
		case 78: goto st78;
		case 84: goto st78;
	}
	goto st0;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	if ( (*p) == 9 )
		goto st58;
	goto st0;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	if ( (*p) == 46 )
		goto st59;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st71;
	goto st0;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	if ( (*p) == 9 )
		goto st60;
	goto st0;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 46: goto st61;
		case 95: goto st69;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st69;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st69;
	} else
		goto st69;
	goto st0;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	if ( (*p) == 9 )
		goto st62;
	goto st0;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	if ( (*p) == 95 )
		goto st63;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st63;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st63;
	} else
		goto st63;
	goto st0;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 9: goto st64;
		case 59: goto st62;
		case 61: goto st67;
		case 95: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st63;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st63;
	} else
		goto st63;
	goto st0;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st65;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st65;
	} else
		goto st65;
	goto st0;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 9: goto st66;
		case 58: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st65;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st65;
	} else
		goto st65;
	goto st0;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st84;
	goto st0;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 9: goto st66;
		case 10: goto st85;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st84;
	goto st0;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	if ( (*p) == 95 )
		goto st49;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st49;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st49;
	} else
		goto st49;
	goto st0;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st68;
	goto st0;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	if ( (*p) == 9 )
		goto st64;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st68;
	goto st0;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 9: goto st62;
		case 44: goto st70;
		case 59: goto st70;
		case 95: goto st69;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st69;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st69;
	} else
		goto st69;
	goto st0;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	if ( (*p) == 95 )
		goto st69;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st69;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st69;
	} else
		goto st69;
	goto st0;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 9: goto st60;
		case 46: goto st72;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st71;
	goto st0;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st73;
	goto st0;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	if ( (*p) == 9 )
		goto st60;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st73;
	goto st0;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	if ( (*p) == 68 )
		goto st75;
	goto st0;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	if ( (*p) == 69 )
		goto st76;
	goto st0;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	if ( (*p) == 76 )
		goto st77;
	goto st0;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	if ( (*p) == 62 )
		goto st57;
	goto st0;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 9: goto st58;
		case 44: goto st79;
		case 65: goto st78;
		case 67: goto st78;
		case 71: goto st78;
		case 78: goto st78;
		case 84: goto st78;
	}
	goto st0;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 65: goto st78;
		case 67: goto st78;
		case 71: goto st78;
		case 78: goto st78;
		case 84: goto st78;
	}
	goto st0;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	if ( (*p) == 9 )
		goto st54;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st80;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st80;
	} else
		goto st80;
	goto st0;
tr88:
#line 38 "vcf.ragel"
	{te = p+1;}
	goto st86;
tr89:
#line 36 "vcf.ragel"
	{te = p+1;{ {cs = stack[--top];goto _again;} }}
	goto st86;
tr91:
#line 32 "vcf.ragel"
	{te = p+1;}
	goto st86;
tr92:
#line 34 "vcf.ragel"
	{te = p+1;}
	goto st86;
tr94:
#line 26 "vcf.ragel"
	{te = p;p--;{ 
            char *fileformat = (char*) calloc (te-ts+1, sizeof(char));
            strncat(fileformat, ts, te-ts);
            set_file_format(fileformat, file);
        }}
	goto st86;
tr103:
#line 24 "vcf.ragel"
	{te = p+1;}
	goto st86;
st86:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof86;
case 86:
#line 1 "NONE"
	{ts = p;}
#line 1066 "vcf_reader.c"
	switch( (*p) ) {
		case 10: goto tr89;
		case 32: goto tr88;
		case 35: goto tr91;
		case 61: goto tr92;
		case 102: goto st88;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto st87;
	} else if ( (*p) >= 9 )
		goto tr88;
	goto st0;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	if ( (*p) == 105 )
		goto st89;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	if ( (*p) == 108 )
		goto st90;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	if ( (*p) == 101 )
		goto st91;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	if ( (*p) == 102 )
		goto st92;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	if ( (*p) == 111 )
		goto st93;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	if ( (*p) == 114 )
		goto st94;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	if ( (*p) == 109 )
		goto st95;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	if ( (*p) == 97 )
		goto st96;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	if ( (*p) == 116 )
		goto tr103;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st87;
	goto tr94;
tr82:
#line 71 "vcf.ragel"
	{{p = ((te))-1;}{
            // Remove quotation marks if necessary
            char *start = ts, *end = te;
            if (*ts == '"') start = ts + 1;
            if (*(te-1) == '"') end = te - 1;
            
            char *field_value = get_token(start, end);
            add_header_entry_value(field_value, status->current_header_entry);
        }}
	goto st97;
tr104:
#line 100 "vcf.ragel"
	{te = p+1;}
	goto st97;
tr105:
#line 93 "vcf.ragel"
	{te = p+1;{
            add_header_entry(status->current_header_entry, file);
            status->current_header_entry = create_header_entry();
            LOG_DEBUG("\n");
            {cs = stack[--top];goto _again;}
        }}
	goto st97;
tr108:
#line 89 "vcf.ragel"
	{te = p+1;{
            LOG_DEBUG(" , ");
        }}
	goto st97;
tr111:
#line 85 "vcf.ragel"
	{te = p+1;{
            LOG_DEBUG(" } ");
        }}
	goto st97;
tr112:
#line 71 "vcf.ragel"
	{te = p;p--;{
            // Remove quotation marks if necessary
            char *start = ts, *end = te;
            if (*ts == '"') start = ts + 1;
            if (*(te-1) == '"') end = te - 1;
            
            char *field_value = get_token(start, end);
            add_header_entry_value(field_value, status->current_header_entry);
        }}
	goto st97;
tr114:
#line 59 "vcf.ragel"
	{te = p;p--;{
            char *field_id = get_token(ts, te-1);
            if (status->current_header_entry->name != NULL)
            {
                add_header_entry_key(field_id, status->current_header_entry);
            } else {
                // Entries like ##reference=some_text_here
                status->current_header_entry = create_header_entry();
                set_header_entry_name(field_id, status->current_header_entry);
            }
        }}
	goto st97;
tr115:
#line 54 "vcf.ragel"
	{te = p+1;{
            char *header_id = get_token(ts, te-2);
            set_header_entry_name(header_id, status->current_header_entry);
        }}
	goto st97;
tr116:
#line 81 "vcf.ragel"
	{te = p+1;{
            LOG_DEBUG(" =< ");
        }}
	goto st97;
st97:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof97;
case 97:
#line 1 "NONE"
	{ts = p;}
#line 1252 "vcf_reader.c"
	switch( (*p) ) {
		case 10: goto tr105;
		case 34: goto tr107;
		case 44: goto tr108;
		case 61: goto st103;
		case 62: goto tr111;
	}
	if ( (*p) < 58 ) {
		if ( (*p) < 32 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr104;
		} else if ( (*p) > 47 ) {
			if ( 48 <= (*p) && (*p) <= 57 )
				goto st101;
		} else
			goto st98;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st101;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st98;
			} else if ( (*p) >= 97 )
				goto st101;
		} else
			goto st98;
	} else
		goto st98;
	goto st0;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	if ( (*p) < 45 ) {
		if ( 32 <= (*p) && (*p) <= 43 )
			goto st98;
	} else if ( (*p) > 61 ) {
		if ( 63 <= (*p) && (*p) <= 126 )
			goto st98;
	} else
		goto st98;
	goto tr112;
tr107:
#line 1 "NONE"
	{te = p+1;}
	goto st99;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
#line 1305 "vcf_reader.c"
	if ( (*p) == 44 )
		goto st81;
	if ( (*p) > 61 ) {
		if ( 63 <= (*p) && (*p) <= 126 )
			goto tr107;
	} else if ( (*p) >= 32 )
		goto tr107;
	goto tr112;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	if ( (*p) == 34 )
		goto tr84;
	if ( (*p) > 61 ) {
		if ( 63 <= (*p) && (*p) <= 126 )
			goto st81;
	} else if ( (*p) >= 32 )
		goto st81;
	goto tr82;
tr84:
#line 1 "NONE"
	{te = p+1;}
	goto st100;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
#line 1334 "vcf_reader.c"
	if ( (*p) == 34 )
		goto tr84;
	if ( (*p) > 61 ) {
		if ( 63 <= (*p) && (*p) <= 126 )
			goto st81;
	} else if ( (*p) >= 32 )
		goto st81;
	goto tr112;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	if ( (*p) == 61 )
		goto st102;
	if ( (*p) < 63 ) {
		if ( (*p) < 45 ) {
			if ( 32 <= (*p) && (*p) <= 43 )
				goto st98;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 60 )
					goto st98;
			} else if ( (*p) >= 48 )
				goto st101;
		} else
			goto st98;
	} else if ( (*p) > 64 ) {
		if ( (*p) < 91 ) {
			if ( 65 <= (*p) && (*p) <= 90 )
				goto st101;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto st98;
			} else if ( (*p) >= 97 )
				goto st101;
		} else
			goto st98;
	} else
		goto st98;
	goto tr112;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	if ( (*p) == 60 )
		goto tr115;
	goto tr114;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	if ( (*p) == 60 )
		goto tr116;
	if ( (*p) < 45 ) {
		if ( 32 <= (*p) && (*p) <= 43 )
			goto st98;
	} else if ( (*p) > 61 ) {
		if ( 63 <= (*p) && (*p) <= 126 )
			goto st98;
	} else
		goto st98;
	goto tr112;
tr117:
#line 125 "vcf.ragel"
	{te = p+1;}
	goto st104;
tr118:
#line 120 "vcf.ragel"
	{te = p+1;{
            LOG_DEBUG("\n");
            {cs = stack[--top];goto _again;} 
        }}
	goto st104;
tr120:
#line 115 "vcf.ragel"
	{te = p;p--;{
            char *sname = get_token(ts, te);
            add_sample_name(sname, file);
        }}
	goto st104;
st104:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof104;
case 104:
#line 1 "NONE"
	{ts = p;}
#line 1424 "vcf_reader.c"
	if ( (*p) == 10 )
		goto tr118;
	if ( (*p) > 13 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st105;
	} else if ( (*p) >= 9 )
		goto tr117;
	goto st0;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st105;
	goto tr120;
tr121:
#line 219 "vcf.ragel"
	{te = p+1;}
	goto st106;
tr122:
#line 198 "vcf.ragel"
	{te = p+1;{
            // If batch is full, add to the list of batches and create a new, empty one
            if (vcf_batch_is_full(status->current_batch))
            {
                list_item_t *item = list_item_new(file->num_records, 1, status->current_batch); 
                list_insert_item(item, batches_list);
                LOG_DEBUG_F("Batch added - %zu records\n", status->current_batch->length);
                status->current_batch = vcf_batch_new(batch_size);
            }

            // If not a blank line, add status->current record to status->current batch
            if (status->current_field != CHROM) {
                add_record_to_vcf_batch(status->current_record, status->current_batch);
                status->num_records++;
            }
            
            status->current_field = CHROM;
            status->num_samples = 0;
            LOG_DEBUG("\n");
        }}
	goto st106;
tr130:
#line 1 "NONE"
	{	switch( act ) {
	case 20:
	{{p = ((te))-1;}
            set_field(ts, te, status);
        }
	break;
	case 22:
	{{p = ((te))-1;}
            set_field(ts, te, status);
        }
	break;
	case 23:
	{{p = ((te))-1;}
            set_field(ts, te, status);
        }
	break;
	case 25:
	{{p = ((te))-1;}
            set_field(ts, te, status);
        }
	break;
	case 27:
	{{p = ((te))-1;}
            set_field(ts, te, status);
        }
	break;
	}
	}
	goto st106;
tr131:
#line 158 "vcf.ragel"
	{te = p;p--;{
            set_field(ts, te, status);
        }}
	goto st106;
tr137:
#line 194 "vcf.ragel"
	{te = p;p--;{
            set_field(ts, te, status);
        }}
	goto st106;
tr139:
#line 182 "vcf.ragel"
	{te = p;p--;{
            set_field(ts, te, status);
        }}
	goto st106;
tr142:
#line 190 "vcf.ragel"
	{te = p;p--;{
            set_field(ts, te, status);
        }}
	goto st106;
tr151:
#line 174 "vcf.ragel"
	{te = p;p--;{
            set_field(ts, te, status);
        }}
	goto st106;
st106:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof106;
case 106:
#line 1 "NONE"
	{ts = p;}
#line 1536 "vcf_reader.c"
	switch( (*p) ) {
		case 10: goto tr122;
		case 32: goto tr121;
		case 46: goto tr124;
		case 60: goto st119;
		case 66: goto st117;
		case 71: goto st123;
		case 78: goto st123;
		case 84: goto st123;
		case 95: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 33 ) {
			if ( 9 <= (*p) && (*p) <= 13 )
				goto tr121;
		} else if ( (*p) > 47 ) {
			if ( (*p) > 57 ) {
				if ( 58 <= (*p) && (*p) <= 64 )
					goto tr123;
			} else if ( (*p) >= 48 )
				goto st108;
		} else
			goto tr123;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st117;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr123;
			} else if ( (*p) >= 97 )
				goto st117;
		} else
			goto tr123;
	} else
		goto st123;
	goto st0;
tr123:
#line 1 "NONE"
	{te = p+1;}
#line 194 "vcf.ragel"
	{act = 27;}
	goto st107;
tr124:
#line 1 "NONE"
	{te = p+1;}
#line 166 "vcf.ragel"
	{act = 20;}
	goto st107;
tr148:
#line 1 "NONE"
	{te = p+1;}
#line 174 "vcf.ragel"
	{act = 22;}
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
#line 1597 "vcf_reader.c"
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr123;
	goto tr130;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	switch( (*p) ) {
		case 44: goto st109;
		case 46: goto tr133;
		case 58: goto st112;
		case 59: goto st114;
		case 61: goto tr136;
		case 95: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 60 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st108;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st117;
	} else
		goto st117;
	goto tr131;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	if ( (*p) == 95 )
		goto st110;
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st110;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st110;
	} else
		goto st110;
	goto tr137;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 44: goto st109;
		case 59: goto st109;
		case 95: goto st110;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st110;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st110;
	} else
		goto st110;
	goto tr139;
tr133:
#line 1 "NONE"
	{te = p+1;}
#line 194 "vcf.ragel"
	{act = 27;}
	goto st111;
tr140:
#line 1 "NONE"
	{te = p+1;}
#line 178 "vcf.ragel"
	{act = 23;}
	goto st111;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
#line 1707 "vcf_reader.c"
	if ( (*p) < 48 ) {
		if ( 33 <= (*p) && (*p) <= 47 )
			goto tr123;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) && (*p) <= 126 )
			goto tr123;
	} else
		goto tr140;
	goto tr130;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st113;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st113;
	} else
		goto st113;
	goto tr137;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	if ( (*p) == 58 )
		goto st112;
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 59 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st113;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st113;
	} else
		goto st113;
	goto tr142;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	if ( (*p) == 95 )
		goto st115;
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st115;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st115;
	} else
		goto st115;
	goto tr137;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 44: goto st109;
		case 59: goto st114;
		case 61: goto tr136;
		case 95: goto st115;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st115;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st115;
	} else
		goto st115;
	goto tr139;
tr136:
#line 1 "NONE"
	{te = p+1;}
#line 194 "vcf.ragel"
	{act = 27;}
	goto st116;
tr144:
#line 1 "NONE"
	{te = p+1;}
#line 186 "vcf.ragel"
	{act = 25;}
	goto st116;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
#line 1843 "vcf_reader.c"
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr144;
	goto tr130;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	switch( (*p) ) {
		case 44: goto st109;
		case 58: goto st112;
		case 59: goto st114;
		case 61: goto tr136;
		case 95: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 60 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st117;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st117;
	} else
		goto st117;
	goto tr131;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 44: goto st109;
		case 59: goto st114;
		case 61: goto tr136;
		case 95: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st118;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto tr123;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) && (*p) <= 126 )
				goto tr123;
		} else
			goto st118;
	} else
		goto st118;
	goto tr131;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	if ( (*p) == 68 )
		goto st120;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr123;
	goto tr137;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	if ( (*p) == 69 )
		goto st121;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr123;
	goto tr137;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	if ( (*p) == 76 )
		goto st122;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr123;
	goto tr137;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	if ( (*p) == 62 )
		goto tr148;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr123;
	goto tr137;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 44: goto st124;
		case 58: goto st112;
		case 59: goto st114;
		case 61: goto tr136;
		case 66: goto st117;
		case 71: goto st123;
		case 78: goto st123;
		case 84: goto st123;
		case 95: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 60 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st117;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st117;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr123;
			} else if ( (*p) >= 97 )
				goto st117;
		} else
			goto tr123;
	} else
		goto st123;
	goto tr131;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 66: goto st110;
		case 71: goto st125;
		case 78: goto st125;
		case 84: goto st125;
		case 95: goto st110;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st110;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st110;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr123;
			} else if ( (*p) >= 97 )
				goto st110;
		} else
			goto tr123;
	} else
		goto st125;
	goto tr137;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 44: goto st124;
		case 59: goto st109;
		case 66: goto st110;
		case 71: goto st125;
		case 78: goto st125;
		case 84: goto st125;
		case 95: goto st110;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 33 <= (*p) && (*p) <= 47 )
				goto tr123;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto tr123;
		} else
			goto st110;
	} else if ( (*p) > 67 ) {
		if ( (*p) < 91 ) {
			if ( 68 <= (*p) && (*p) <= 90 )
				goto st110;
		} else if ( (*p) > 96 ) {
			if ( (*p) > 122 ) {
				if ( 123 <= (*p) && (*p) <= 126 )
					goto tr123;
			} else if ( (*p) >= 97 )
				goto st110;
		} else
			goto tr123;
	} else
		goto st125;
	goto tr151;
	}
	_test_eof82: cs = 82; goto _test_eof; 
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
	_test_eof83: cs = 83; goto _test_eof; 
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
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
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
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
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
	_test_eof81: cs = 81; goto _test_eof; 
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

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 87: goto tr94;
	case 88: goto tr94;
	case 89: goto tr94;
	case 90: goto tr94;
	case 91: goto tr94;
	case 92: goto tr94;
	case 93: goto tr94;
	case 94: goto tr94;
	case 95: goto tr94;
	case 96: goto tr94;
	case 98: goto tr112;
	case 99: goto tr112;
	case 81: goto tr82;
	case 100: goto tr112;
	case 101: goto tr112;
	case 102: goto tr114;
	case 103: goto tr112;
	case 105: goto tr120;
	case 107: goto tr130;
	case 108: goto tr131;
	case 109: goto tr137;
	case 110: goto tr139;
	case 111: goto tr130;
	case 112: goto tr137;
	case 113: goto tr142;
	case 114: goto tr137;
	case 115: goto tr139;
	case 116: goto tr130;
	case 117: goto tr131;
	case 118: goto tr131;
	case 119: goto tr137;
	case 120: goto tr137;
	case 121: goto tr137;
	case 122: goto tr137;
	case 123: goto tr131;
	case 124: goto tr137;
	case 125: goto tr151;
	}
	}

	_out: {}
	}

#line 302 "vcf.ragel"

    
    if (status->self_contained && !vcf_batch_is_empty(status->current_batch)) {
        list_item_t *item = list_item_new(file->num_records, 1, status->current_batch);
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (self-contained)\n", status->current_batch->length);
    }

    return cs;
}


int vcf_ragel_read(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples) {
    int cs;
    char *p, *pe;

    vcf_reader_status *status = new_vcf_reader_status(batch_size, read_samples, 0);
    
    if (mmap_vcf) {
        LOG_DEBUG("Using mmap for file loading\n");
        p = file->data;
        pe = p + file->data_len;
        cs = execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
    } else {
        LOG_DEBUG("Using file-IO functions for file loading\n");
        size_t max_len = 256;
        __ssize_t line_len = 0;

        char *data = (char*) malloc (max_len * sizeof(char));
        char *line = NULL;
        char *aux;
        
        int eof_found = 0;

        // Read text of a batch and call ragel parser in a loop
        while (!eof_found) {
            memset(data, 0, max_len * sizeof(char));

            for (int i = 0; i < batch_size && !eof_found; i++) {
                line_len = getline(&line, &line_len, file->fd);
                if (line_len != -1) {
    //                LOG_DEBUG_F("Line (len %zu): %s", line_len, line);
                    // Line too long to be stored in data, realloc
                    if (file->data_len + line_len + 1 > max_len) {
                        aux = realloc(data, max_len + line_len * 20);
                        if (aux) {
                            data = aux;
                            max_len += line_len * 20;
                        } else {
                            LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
                        }
                    }
                    // Concat previous data with new line
                    strncat(data, line, line_len);
                    file->data_len += line_len;
                } else {
                    eof_found = 1;
                }
            }

            p = data;
            pe = p + file->data_len;
            cs = execute_vcf_ragel_machine(p, pe, batches_list, batch_size, file, status);
            file->data_len = 0;
        }
        
        if (line != NULL) { free(line); }
        if (data != NULL) { free(data); }
    }

    // Insert the last batch
    if (!vcf_batch_is_empty(status->current_batch))
    {
        list_item_t *item = list_item_new(file->num_records, 1, status->current_batch);
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", status->current_batch->length);
    }

    if ( cs < 
#line 2311 "vcf_reader.c"
82
#line 380 "vcf.ragel"
 ) 
    {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 2318 "vcf_reader.c"
82
#line 383 "vcf.ragel"
);
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", file->num_samples);

    // Free status->current_xxx pointers if not needed in another module
    vcf_header_entry_free(status->current_header_entry);   

    return cs < 
#line 2331 "vcf_reader.c"
82
#line 392 "vcf.ragel"
;
}

int vcf_light_read(list_t *batches_list, size_t batch_size, vcf_file_t *file) {
//     if (mmap_vcf) {
//         LOG_DEBUG("Using mmap for file loading\n");
//         TODO Use functions for getting the batch-size-th index of linebreak
//     } else {
    LOG_DEBUG("Using file-IO functions for file loading\n");
    size_t max_len = 256;
    __ssize_t line_len = 0;
    file->data_len = 0;

    char *data = NULL;
    char *line = NULL;
    char *aux;
    
    int eof_found = 0;

    // Read text of a batch and call ragel parser in a loop
    while (!eof_found) {
        char *data = (char*) calloc (max_len, sizeof(char));
//             memset(data, 0, max_len * sizeof(char));

        for (int i = 0; i < batch_size && !eof_found; i++) {
            line_len = getline(&line, &line_len, file->fd);
            if (line_len != -1) {
                LOG_DEBUG_F("#%d Line (len %zu): %s", i, line_len, line);
                // Line too long to be stored in data, realloc
                if (file->data_len + line_len + 1 > max_len) {
                    aux = realloc(data, max_len + line_len * 20);
                    if (aux) {
                        data = aux;
                        max_len += line_len * 20;
                    } else {
                        LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
                    }
                }
                // Concat previous data with new line
                strncat(data, line, line_len);
                file->data_len += line_len;
            } else {
                eof_found = 1;
            }
        }

        file->data_len = 0;

        list_item_t *item = list_item_new(file->num_records, 1, data);
        list_insert_item(item, batches_list);
//             printf("Text batch inserted = '%s'\n", data);

    }
    
    if (line != NULL) { free(line); }

    return 0;
}

vcf_reader_status *new_vcf_reader_status(size_t batch_size, int store_samples, int self_contained) {
    vcf_reader_status *status = (vcf_reader_status *) malloc (sizeof(vcf_reader_status));
    status->current_record = NULL;
    status->current_header_entry = create_header_entry();
    status->current_batch = vcf_batch_new(batch_size);
    status->current_field = CHROM;

    status->num_samples = 0;
    status->num_records = 0;

    status->store_samples = store_samples;
    status->self_contained = self_contained;

    return status;
}
