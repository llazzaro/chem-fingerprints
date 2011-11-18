#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>

#include "chemfp.h"
#include "chemfp_internal.h"
#include "popcount.h"

#include "cpuid.h"

//#define REGULAR

chemfp_alignment_type _chemfp_alignments[] = {
  {"align1", 1, 1, NULL},
  {"align4", 4, 4, NULL},
  {"align8-small", 8, 8, NULL},
  {"align8-large", 8, 96, NULL},
};

static int
has_popcnt_instruction(void) {
  return (get_cpuid_flags() & bit_POPCNT);
}

static chemfp_method_type compile_time_methods[] = {
  {0, "LUT8-1", 1, 1, 0, NULL,
   chemfp_byte_popcount,
    chemfp_byte_intersect_popcount},

  {0, "LUT8-4", 4, 4, 1, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut8_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut8_4},

  {0, "LUT16-4", 4, 4, 1, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut16_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut16_4},

  {0, "Lauradoux", 8, 96, 1, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lauradoux,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lauradoux},

  {0, "popcnt", 8, 8, 1,
   has_popcnt_instruction,
   (chemfp_popcount_f) popcount_POPCNT,
   (chemfp_intersect_popcount_f) intersect_popcount_POPCNT},

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

static inline void
set_alignment_methods(void) {
  int alignment;

  /* Make sure we haven't already initialized the alignments */
  if (_chemfp_alignments[0].method_p != NULL) {
    return;
  }

  /* Figure out which methods are available for this hardware */
  detect_methods();

  /* Initialize to something which is always valid */
  for (alignment=0; alignment < sizeof(_chemfp_alignments) / sizeof(chemfp_alignment_type);
       alignment++) {
    _chemfp_alignments[alignment].method_p = detected_methods[0];
  }

  /* Determine the fastest method for the different alignment categories */
  /* (I chose "1000" because the times end up around 100 microseconds. */
  /*  Any lower and the fastest times might be too small for high-end machines) */
  chemfp_select_fastest_method(CHEMFP_ALIGN4, 1000);
  chemfp_select_fastest_method(CHEMFP_ALIGN8_SMALL, 1000);
  chemfp_select_fastest_method(CHEMFP_ALIGN8_LARGE, 1000);
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


/*********** Automatically select the fastest method ***********/

static const int buffersize = 2048/8;

static unsigned long
get_usecs(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec*1000000+tv.tv_usec;
}

static unsigned long 
timeit(chemfp_popcount_f popcount, int size, unsigned char *buffer, int repeat) {
  unsigned long t1, t2;
  int i;
  t1 = get_usecs();
  for (i=0; i<repeat; i++) {
    popcount(size, buffer);
  }
  t2 = get_usecs();
  return t2-t1;
}



int
chemfp_select_fastest_method(int alignment, int repeat) {
  unsigned char *buffer;
  int method, best_method=-1, old_method;
  int i, probe_size;
  unsigned long dt;
  unsigned long best_time=0;
  chemfp_method_type *method_p=NULL;

  old_method = chemfp_get_alignment_method(alignment);
  if (old_method < 0) {
    return old_method;
  }

  buffer = (unsigned char *) malloc(buffersize);
  for (i=0; i<buffersize; i++) {
    buffer[i] = i % 256;
  }

  
  if (alignment == CHEMFP_ALIGN8_SMALL) {
    probe_size = 64; /* 512 bits */
  } else {
    probe_size = 2048/8;
  }

  for (method=0; method<chemfp_get_num_methods(); method++) {

    /* See if I can use this method */
    if (chemfp_set_alignment_method(alignment, method) < 0) {
      continue;
    }
    method_p = _chemfp_alignments[alignment].method_p;
    if (method_p->alignment < _chemfp_alignments[alignment].alignment &&
	!method_p->okay_for_larger_alignments) {
      continue;
    }

    /* Time the performance */
    dt = timeit(method_p->popcount,
		probe_size, buffer, repeat);
		
    if (best_method == -1 || dt < best_time) {
      best_method = method;
      best_time = dt;
    }
  }
  if (best_method == -1) {
    /* Shouldn't happen, but I want to be on the safe side. */
    best_method = old_method;
  }
  chemfp_set_alignment_method(alignment, best_method);

  free(buffer);

  return best_method;
}
