#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"
#include "popcount.h"

//#define REGULAR

chemfp_alignment_type _chemfp_alignments[] = {
  {"align1", 1, 1, NULL},
  {"align4", 4, 4, NULL},
  {"align8-small", 8, 8, NULL},
  {"align8-large", 8, 96, NULL},
};


static chemfp_method_type compile_time_methods[] = {
  {0, "LUT8-1", 1, 1, NULL,
   chemfp_byte_popcount,
    chemfp_byte_intersect_popcount},

  {0, "LUT8-4", 4, 4, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut8_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut8_4},

  {0, "LUT16-4", 4, 4, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut16_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut16_4},

  {0, "Lauradoux", 8, 96, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lauradoux,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lauradoux},

#if defined(HAS_POPCOUNT_INTRINSIC)
  {0, "intrinsic32", 4, 4, /* To implement */
   _chemfp_has_popcnt_instruction,
   (chemfp_popcount_f) _chemfp_popcount_intrinsic32,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_intrinsic32},

  {0, "intrinsic64", 8, 8, /* To implement */
   _chemfp_has_popcnt_instruction,
   (chemfp_popcount_f) _chemfp_popcount_intrinsic64,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_intrinsic64},
#endif

};


/* These are the methods which are actually available at run-time */

static chemfp_method_type *
detected_methods[sizeof(compile_time_methods)/sizeof(chemfp_method_type)];

static int num_methods = 0;

static void
detect_methods(void) {
  int i, j=0;
  if (num_methods != 0) {
    return;
  }
  for (i=0; i<sizeof(compile_time_methods)/sizeof(chemfp_method_type); i++) {
    if ((compile_time_methods[i].check == NULL) ||
	(compile_time_methods[i].check())) {
      compile_time_methods[i].detected_index = j;
      detected_methods[j++] = &compile_time_methods[i];
    }
  }
  num_methods = j;
}



int
chemfp_get_num_methods(void) {
  detect_methods();
  return num_methods;
}

const char *
chemfp_get_method_name(int method) {
  if (method < 0 || method >= chemfp_get_num_methods()) {
    return NULL;
  }
  return detected_methods[method]->name;
}

static chemfp_method_type *
_get_method(const char **names, int num_bytes) {
  int probe, i;
  int num_names = num_bytes / sizeof(const char *);
  for (probe=0; probe<num_names; probe++) {
    for (i=0; i<num_methods; i++) {
      if (!strcmp(names[probe], detected_methods[i]->name)) {
	printf("Assign %s\n", detected_methods[i]->name);
	return detected_methods[i];
      }
    }
  }
  printf("Default\n");
  return &compile_time_methods[0];
}

	    
static inline void
set_alignment_methods(void) {
  if (_chemfp_alignments[0].method_p != NULL) {
    return;
  }
  printf("detect methods\n");
  detect_methods();
  printf("detect methods is done\n");
  const char *align1[] = {"LUT8-1"};
  const char *align4[] = {"LUT16-4", "LUT8-4"};
  const char *align8_small[] = {"intrinsic64", "LUT16-4", "LUT8-4"};
  const char *align8_large[] = {"intrinsic64", "Lauradoux", "LUT16-4", "LUT8-4"};
  _chemfp_alignments[0].method_p = _get_method(align1, sizeof(align1));
  _chemfp_alignments[1].method_p = _get_method(align4, sizeof(align4));
  _chemfp_alignments[2].method_p = _get_method(align8_small, sizeof(align8_small));
  _chemfp_alignments[3].method_p = _get_method(align8_large, sizeof(align8_large));
}


int
chemfp_get_num_alignments(void) {
  set_alignment_methods();
  return sizeof(_chemfp_alignments) / sizeof(chemfp_alignment_type);
}

const char *
chemfp_get_alignment_name(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return NULL;
  }
  return _chemfp_alignments[alignment].name;
}

int 
chemfp_get_alignment_method(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return CHEMFP_BAD_ARG;
  }
  return _chemfp_alignments[alignment].method_p->detected_index;
}

int
chemfp_set_alignment_method(int alignment, int method) {
  /* Make sure it's an available alignment and method */
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return CHEMFP_BAD_ARG;
  }
  if (method < 0 || method >= chemfp_get_num_methods()) {
    return CHEMFP_BAD_ARG;
  }
  /* Make sure the alignment and sizes are good enough */
  if (detected_methods[method]->alignment > _chemfp_alignments[alignment].alignment) {
    return CHEMFP_METHOD_MISMATCH;
  }
  if (detected_methods[method]->min_size > _chemfp_alignments[alignment].min_size) {
    return CHEMFP_METHOD_MISMATCH;
  }
  _chemfp_alignments[alignment].method_p = detected_methods[method];
  return CHEMFP_OK;
}

  

/**************************************/

/* chemfp stores fingerprints as Python strings */
/* (This may change in the future; memmap, perhaps?) */
/* The Python payload is 4 byte aligned but not 8 byte aligned. */

chemfp_popcount_f
chemfp_select_popcount(int num_bits,
		       int storage_len, const unsigned char *arena) {

  int num_bytes = (num_bits+7)/8;

#ifdef REGULAR
  return chemfp_byte_popcount;
#endif
  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_alignment_methods();

  if (num_bytes <= 1) {
    /* Really? */
    return _chemfp_alignments[CHEMFP_ALIGN1].method_p->popcount;
  }
  if (ALIGNMENT(arena, 8) == 0 &&
      storage_len % 8 == 0) {
    if (num_bytes >= 96) {
      return _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p->popcount;
    } else {
      return _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->popcount;
    }
  }

  if (ALIGNMENT(arena, 4) &&
      storage_len % 4 == 0) {
    return _chemfp_alignments[CHEMFP_ALIGN4].method_p->popcount;
  }

  return _chemfp_alignments[CHEMFP_ALIGN1].method_p->popcount;
}



chemfp_intersect_popcount_f
chemfp_select_intersect_popcount(int num_bits,
				 int storage_len1, const unsigned char *arena1,
				 int storage_len2, const unsigned char *arena2) {

  int storage_len = (storage_len1 < storage_len2) ? storage_len1 : storage_len2;

  int num_bytes = (num_bits+7)/8;

#ifdef REGULAR
  return chemfp_byte_intersect_popcount;
#endif

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_alignment_methods();
  
  if (num_bytes <= 1) {
    return _chemfp_alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
  }

  /* Check for 8 byte alignment */

  if (ALIGNMENT(arena1, 8) == 0 &&
      ALIGNMENT(arena2, 8) == 0 &&
      storage_len1 % 8 == 0 &&
      storage_len2 % 8 == 0) {

    if (num_bytes >= 96) {
      return _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p->intersect_popcount;
    } else {
      return _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->intersect_popcount;
    }
  }

  /* Check for 4 byte alignment */

  if (ALIGNMENT(arena1, 4) == 0 &&
      ALIGNMENT(arena2, 4) == 0 &&
      storage_len1 % 4 == 0 &&
      storage_len2 % 4 == 0) {
    return _chemfp_alignments[CHEMFP_ALIGN4].method_p->intersect_popcount;
  }

  /* At least we're one byte aligned */
  return _chemfp_alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
}
