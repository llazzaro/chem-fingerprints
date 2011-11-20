#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>

#include "chemfp.h"
#include "chemfp_internal.h"
#include "popcount.h"

#include "cpuid.h"

static unsigned long 
timeit(chemfp_popcount_f popcount, int size, int repeat);

/* These are the alignment categories which I support */

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

/* These are in the same order as an enum in popcount.h */
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

  {0, "popcnt", 8, 8,
   has_popcnt_instruction,
   (chemfp_popcount_f) popcount_POPCNT,
   (chemfp_intersect_popcount_f) intersect_popcount_POPCNT},

};


/* These are the methods which are actually available at run-time */
/* This list is used for the public API */

static chemfp_method_type *
detected_methods[sizeof(compile_time_methods)/sizeof(chemfp_method_type)];

static int num_methods = 0;

static void
detect_methods(void) {
  int i, j=0;
  if (num_methods != 0) {
    return;
  }
  /* Go through all of the compile-time methods and see if it's available */
  for (i=0; i<sizeof(compile_time_methods)/sizeof(chemfp_method_type); i++) {
    if ((compile_time_methods[i].check == NULL) ||
	(compile_time_methods[i].check())) {

      /* Add it to the list of detected methods, and tell it its index position */
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
set_default_alignment_methods(void) {
  int lut_method, large_method;
  unsigned long first_time, lut8_time, lut16_time, lut_time, lauradoux_time;

  /* Make sure we haven't already initialized the alignments */
  if (_chemfp_alignments[0].method_p != NULL) {
    return;
  }

  /* Figure out which methods are available for this hardware */
  detect_methods();

  /* This is the only possibility for 1-byte aligned */
  _chemfp_alignments[CHEMFP_ALIGN1].method_p = &compile_time_methods[CHEMFP_LUT8_1];

  /* Now do some timing measurements and figure out which method is
     likely the fastest for this hardware. It's a bit tricky; consider
     what happens if a timeslice boundary happens while doing a
     test. I mostly fix that by doing the timing twice and using
     the fastest time.

     I could require everyone call chemfp_select_fastest_method,
     but this should be good enough for almost everyone. */


  /* For 4-byte aligned we use a LUT. */
  /* TODO: implement a POPCNT instruction-based method for 4-byte aligned code */
  /* (You really should use an 8-byte aligned arena in this case, so not a priority */

  /* On older hardware the LUT16 can be slower than the LUT8 */
  first_time = timeit(compile_time_methods[CHEMFP_LUT8_4].popcount, 128, 700);
  lut8_time = timeit(compile_time_methods[CHEMFP_LUT8_4].popcount, 128, 700);
  if (first_time < lut8_time) {
    lut8_time = first_time;
  }

  first_time = timeit(compile_time_methods[CHEMFP_LUT16_4].popcount, 128, 700);
  lut16_time = timeit(compile_time_methods[CHEMFP_LUT16_4].popcount, 128, 700);
  if (first_time < lut16_time) {
    lut16_time = first_time;
  }

  /* Which one is faster? */
  if (lut8_time < lut16_time) {
    lut_method = CHEMFP_LUT8_4;
    lut_time = lut8_time;
  } else {
    lut_method = CHEMFP_LUT16_4;
    lut_time = lut16_time;
  }
  _chemfp_alignments[CHEMFP_ALIGN4].method_p = &compile_time_methods[lut_method];


  /* For 8-byte aligned code we always want to use the POPCNT instruction if it exists */
  if (has_popcnt_instruction()) {
    _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p = 
      _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p = &compile_time_methods[CHEMFP_POPCNT];
  } else {

    /* No POPCNT? Then it's a LUT for the small case, and perhaps Lauradoux */
    /* for the large case */

    _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p = &compile_time_methods[lut_method];

    first_time = timeit(compile_time_methods[CHEMFP_LAURADOUX].popcount, 128, 700);
    lauradoux_time = timeit(compile_time_methods[CHEMFP_LAURADOUX].popcount, 128, 700);
    if (first_time < lauradoux_time) {
      lauradoux_time = first_time;
    }

    if (lauradoux_time < lut_time) {
      large_method = CHEMFP_LAURADOUX;
    } else {
      large_method = lut_method;
    }
    _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p = &compile_time_methods[large_method];
  }
}


int
chemfp_get_num_alignments(void) {
  set_default_alignment_methods();
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

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_default_alignment_methods();

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

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_default_alignment_methods();
  
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

static unsigned long
get_usecs(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec*1000000+tv.tv_usec;
}


/* Space for 2048 bits == 32 * 64 bit words */
/* Use uint64_t so it's 64-bit/8 byte aligned */
/* The contents are randomly generated. */
static uint64_t popcount_buffer[] = {
  0x9b649615d1a50133ull,
  0xf3b8dada0e8b43deull,
  0x0197e207e4b9af2bull,
  0x68a2ecc4053b1305ull,
  0x93d933ac2f41e28full,
  0xb460859e01b6f925ull,
  0xc2c1a9eacc9e4999ull,
  0xdc5237f8200aec07ull,
  0x9e3bbe45d6e67641ull,
  0xa49bed7d060407d4ull,
  0xcca5f2913af53c5bull,
  0xfdd53575aab7c21aull,
  0x76b82d57bfa5c9ddull,
  0x0d2a87ba7f2439edull,
  0x9ec6a4ee2a6999d4ull,
  0xb9ae55f1f402ac97ull,
  0x08bbc6d1719a56bdull,
  0x969e5ef023c9ed23ull,
  0x6b7f08af661a9db6ull,
  0xad394da52bbbe18dull,
  0xdf9c3e28aae1c460ull,
  0xcf82e77d4f02f1efull,
  0x1fb88cdb648008ecull,
  0xc7a2ab7ecb8f84f5ull,
  0xbf8ef6833f18d407ull,
  0xb9c7eafdb4653fa2ull,
  0x90114b93b87a8a1dull,
  0x6e572c9e42e5061cull,
  0xb694ec549eeabc20ull,
  0xb362909621b9a2c8ull,
  0xcadab7b921d3cd0aull,
  0xd27f7aef7e2a0c6full,
};

static unsigned long 
timeit(chemfp_popcount_f popcount, int size, int repeat) {
  unsigned long t1, t2;
  int i;
  if (size > sizeof(popcount_buffer)) {
    size = sizeof(popcount_buffer);
  }
  t1 = get_usecs();
  for (i=0; i<repeat; i++) {
    popcount(size, (unsigned char*) popcount_buffer);
  }
  t2 = get_usecs();
  return t2-t1;
}



int
chemfp_select_fastest_method(int alignment, int repeat) {
  int method, best_method=-1, old_method;
  int probe_size;
  unsigned long dt;
  unsigned long best_time=0;
  chemfp_method_type *method_p=NULL;

  old_method = chemfp_get_alignment_method(alignment);
  if (old_method < 0) {
    return old_method;
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

    /* Time the performance */
    dt = timeit(method_p->popcount, probe_size, repeat);
		
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

  return best_method;
}
