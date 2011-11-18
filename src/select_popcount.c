#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"

//#define REGULAR

#include "lut16.c"


/* This is called with the number of bytes in the fingerprint, */
/* which is not necessarily a multiple of 4. However, the select */
/* function promises that there will be enough extra data. */

static inline int popcount_lut8(int n, uint32_t *fp) {
  int cnt=0;
  unsigned int i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    i = *fp;
    cnt += lut[i&255];
    cnt += lut[i>>8&255];
    cnt += lut[i>>16&255];
    cnt += lut[i>>24];
    fp++;
  } while(--n);
  return cnt;
}

static inline int intersect_popcount_lut8(int n, uint32_t *fp1, uint32_t *fp2) {
  int cnt=0;
  unsigned int i;
  n = (n+3) / 4;
  do {
    i = *fp1 & *fp2;
    cnt += lut[i&255];
    cnt += lut[i>>8&255];
    cnt += lut[i>>16&255];
    cnt += lut[i>>24];  /* APD: removed the unneeded &255 */
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}

static inline int popcount_lut16(int n, uint32_t *fp) {
  int cnt=0;
  unsigned int i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    i = *fp;
    cnt += lut[i&65535];
    cnt += lut[i>>16];
    fp++;
  } while(--n);
  return cnt;
}

static inline int intersect_popcount_lut16(int n, uint32_t *fp1, uint32_t *fp2) {
  int cnt=0;
  unsigned int i;
  n = (n+3) / 4;
  do {
    i = *fp1 & *fp2;
    cnt += lut[i&65535];
    cnt += lut[i>>16];
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}


/**************************************/

#if defined(HAS_POPCOUNT_INTRINSIC)

static inline int popcount_intrinsic32(int n, uint32_t *fp) {
  int cnt=0;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    cnt += POPCNT32(*fp++);
  } while(--n);
  return cnt;
}

static inline int intersect_popcount_intrinsic32(int n, uint32_t *fp1, uint32_t *fp2) {
  int cnt=0;
  uint32_t i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    i = (*fp1 & *fp2);
    cnt += POPCNT32(i);
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}


static inline int popcount_intrinsic64(int n, uint64_t *fp) {
  int cnt=0;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+7) / 8;
  do {
    cnt += POPCNT64(*fp++);
  } while(--n);
  return cnt;
}

static inline int intersect_popcount_intrinsic64(int n, uint64_t *fp1, uint64_t *fp2) {
  int cnt=0;
  uint64_t i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+7) / 8;
  do {
    i = (*fp1 & *fp2);
    cnt += POPCNT64(i);
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}

static int
has_popcnt(void) {
  printf("Do I have popcount\n");
  return 1;
}
#endif


/**************************************/
//
// popcount_lauradoux.c, v 1.1
// Created by Kim Walisch on 11/10/11.
//
// Changes:
//
// The main outer loop of the two functions now processes
// 12*8 = 96 bytes per iteration (previously 240 bytes).
// This makes the functions more efficient for small fingerprints
// e.g. 881 bits.

#include <assert.h>

/**
 * Count the number of 1 bits (population count) in an array using
 * 64-bit tree merging. This implementation uses only 8 operations per
 * 8 bytes on 64-bit CPUs.
 * The algorithm is due to Cédric Lauradoux, it is described and
 * benchmarked against other bit population count solutions (lookup
 * tables, bit-slicing) in his paper:
 * http://perso.citi.insa-lyon.fr/claurado/ham/overview.pdf
 * http://perso.citi.insa-lyon.fr/claurado/hamming.html
 */
static int
popcount_lauradoux(int size, const uint64_t *fp) {
  assert(fp != NULL);
  assert(size <= UINT32_MAX / (8 * sizeof(uint64_t)));
  
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
  const uint64_t h01 = UINT64_C(0x0101010101010101);
  
  uint64_t count1, count2, half1, half2, acc;
  uint64_t x;
  int i, j;    
  size = (size + 7) / 8;
  int limit = size - size % 12;
  int bit_count = 0;
  
  // 64-bit tree merging
  for (i = 0; i < limit; i += 12, fp += 12) {
    acc = 0;
    for (j = 0; j < 12; j += 3) {
      count1  =  fp[j + 0];
      count2  =  fp[j + 1];
      half1   =  fp[j + 2];
      half2   =  fp[j + 2];
      half1  &=  m1;
      half2   = (half2  >> 1) & m1;
      count1 -= (count1 >> 1) & m1;
      count2 -= (count2 >> 1) & m1;
      count1 +=  half1;
      count2 +=  half2;
      count1  = (count1 & m2) + ((count1 >> 2) & m2);
      count1 += (count2 & m2) + ((count2 >> 2) & m2);
      acc    += (count1 & m4) + ((count1 >> 4) & m4);
    }
    acc = (acc & m8) + ((acc >>  8)  & m8);
    acc = (acc       +  (acc >> 16)) & m16;
    acc =  acc       +  (acc >> 32);
    bit_count += (int) acc;
  }
  
  // count the bits of the remaining bytes (MAX 88) using 
  // "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
  // the code uses wikipedia's 64-bit popcount_3() implementation:
  // http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
  for (i = 0; i < size - limit; i++) {
    x = fp[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}

static int
intersect_popcount_lauradoux( int size, const uint64_t *fp1, const uint64_t *fp2) {
  assert(fp1 != NULL && fp2 != NULL);
  assert(size <= UINT32_MAX / (8 * sizeof(uint64_t)));
  
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
  const uint64_t h01 = UINT64_C(0x0101010101010101);
  
  uint64_t count1, count2, half1, half2, acc;
  uint64_t x;
  int i, j;
  /* Adjust the input size because the size isn't necessarily a multiple of 8 */
  /* even though the storage area will be a multiple of 8. */
  size = (size + 7) / 8;
  int limit = size - size % 12;
  int bit_count = 0;
  
  // 64-bit tree merging
  for (i = 0; i < limit; i += 12, fp1 += 12, fp2 += 12) {
    acc = 0;
    for (j = 0; j < 12; j += 3) {
      count1  =  fp1[j + 0] & fp2[j + 0];
      count2  =  fp1[j + 1] & fp2[j + 1];
      half1   =  fp1[j + 2] & fp2[j + 2];
      half2   =  fp1[j + 2] & fp2[j + 2];
      half1  &=  m1;
      half2   = (half2  >> 1) & m1;
      count1 -= (count1 >> 1) & m1;
      count2 -= (count2 >> 1) & m1;
      count1 +=  half1;
      count2 +=  half2;
      count1  = (count1 & m2) + ((count1 >> 2) & m2);
      count1 += (count2 & m2) + ((count2 >> 2) & m2);
      acc    += (count1 & m4) + ((count1 >> 4) & m4);
    }
    acc = (acc & m8) + ((acc >>  8)  & m8);
    acc = (acc       +  (acc >> 16)) & m16;
    acc =  acc       +  (acc >> 32);
    bit_count += (int) acc;
  }
  
  // intersect count the bits of the remaining bytes (MAX 88) using 
  // "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
  // the code uses wikipedia's 64-bit popcount_3() implementation:
  // http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
  for (i = 0; i < size - limit; i++) {
    x = fp1[i] & fp2[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}


/************************/

typedef int (*method_check_f)(void);


typedef struct {
  int detected_index;
  const char *name;
  int alignment;
  int min_size;
  method_check_f check;
  chemfp_popcount_f popcount;
  chemfp_intersect_popcount_f intersect_popcount;
} method_type;

static method_type available_methods[] = {
  {0, "LUT8-1", 1, 1,
   NULL, chemfp_byte_popcount, chemfp_byte_intersect_popcount},

  {0, "LUT8-4", 4, 4,
   NULL, popcount_lut8, intersect_popcount_lut8},

  {0, "LUT16-4", 4, 4,  /* To implement */
   NULL, popcount_lut16, intersect_popcount_lut16},

#if defined(HAS_POPCOUNT_INTRINSIC)
  {0, "intrinsic32", 4, 4, /* To implement */
   has_popcnt, popcount_intrinsic32, intersect_popcount_intrinsic32},

  {0, "intrinsic64", 8, 8, /* To implement */
   has_popcnt, popcount_intrinsic64, intersect_popcount_intrinsic64},
#endif

  {0, "Lauradoux", 8, 96,
   NULL, popcount_lauradoux, intersect_popcount_lauradoux},
};

static method_type *methods[sizeof(available_methods)/sizeof(method_type)];
static int num_methods = 0;

static void
detect_methods(void) {
  int i, j=0;
  if (num_methods != 0) {
    return;
  }
  for (i=0; i<sizeof(available_methods)/sizeof(method_type); i++) {
    if ((available_methods[i].check == NULL) ||
	(available_methods[i].check())) {
      available_methods[i].detected_index = j;
      methods[j++] = &available_methods[i];
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
  return methods[method]->name;
}


typedef struct {
  const char *name;
  int alignment;
  int min_size;
  method_type *method_p;
} alignment_type;

static alignment_type alignments[] = {
  {"align1", 1, 1, NULL},
  {"align4", 4, 4, NULL},
  {"align8-small", 8, 8, NULL},
  {"align8-large", 8, 96, NULL},
};
enum {
  CHEMFP_ALIGN1=0,
  CHEMFP_ALIGN4,
  CHEMFP_ALIGN8_SMALL,
  CHEMFP_ALIGN8_LARGE,
};


static method_type *
_get_method(const char **names, int num_bytes) {
  int probe, i;
  int num_names = num_bytes / sizeof(const char *);
  for (probe=0; probe<num_names; probe++) {
    for (i=0; i<num_methods; i++) {
      if (!strcmp(names[probe], methods[i]->name)) {
	printf("Assign %s\n", methods[i]->name);
	return methods[i];
      }
    }
  }
  printf("Default\n");
  return &available_methods[0];
}

	    
static inline void
set_alignment_methods(void) {
  if (alignments[0].method_p != NULL) {
    return;
  }
  printf("detect methods\n");
  detect_methods();
  printf("detect methods is done\n");
  const char *align1[] = {"LUT8-1"};
  const char *align4[] = {"LUT16-4", "LUT8-4"};
  const char *align8_small[] = {"intrinsic64", "LUT16-4", "LUT8-4"};
  const char *align8_large[] = {"intrinsic64", "Lauradoux", "LUT16-4", "LUT8-4"};
  alignments[0].method_p = _get_method(align1, sizeof(align1));
  alignments[1].method_p = _get_method(align4, sizeof(align4));
  alignments[2].method_p = _get_method(align8_small, sizeof(align8_small));
  alignments[3].method_p = _get_method(align8_large, sizeof(align8_large));
}


int
chemfp_get_num_alignments(void) {
  set_alignment_methods();
  return sizeof(alignments) / sizeof(alignment_type);
}

const char *
chemfp_get_alignment_name(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return NULL;
  }
  return alignments[alignment].name;
}

int 
chemfp_get_alignment_method(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return -1;
  }
  return alignments[alignment].method_p->detected_index;
}

int
chemfp_set_alignment_method(int alignment, int method) {
  /* Make sure it's an available alignment and method */
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return -1;
  }
  if (method < 0 || method >= chemfp_get_num_methods()) {
    return -1;
  }
  /* Make sure the alignment and sizes are good enough */
  if (methods[method]->alignment > alignments[alignment].alignment) {
    return -1;
  }
  if (methods[method]->min_size > alignments[alignment].min_size) {
    return -1;
  }
  alignments[alignment].method_p = methods[method];
  return 0;
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
    return alignments[CHEMFP_ALIGN1].method_p->popcount;
  }
  if (ALIGNMENT(arena, 8) == 0 &&
      storage_len % 8 == 0) {
    if (num_bytes >= 96) {
      return alignments[CHEMFP_ALIGN8_LARGE].method_p->popcount;
    } else {
      return alignments[CHEMFP_ALIGN8_SMALL].method_p->popcount;
    }
  }

  if (ALIGNMENT(arena, 4) &&
      storage_len % 4 == 0) {
    return alignments[CHEMFP_ALIGN4].method_p->popcount;
  }

  if (ALIGNMENT(arena, 2) &&
      storage_len % 2 == 0) {
    return alignments[CHEMFP_ALIGN1].method_p->popcount;
  }

  return alignments[CHEMFP_ALIGN1].method_p->popcount;
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
    return alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
  }

  /* Check for 8 byte alignment */

  if (ALIGNMENT(arena1, 8) == 0 &&
      ALIGNMENT(arena2, 8) == 0 &&
      storage_len1 % 8 == 0 &&
      storage_len2 % 8 == 0) {

    if (num_bytes >= 96) {
      return alignments[CHEMFP_ALIGN8_LARGE].method_p->intersect_popcount;
    } else {
      return alignments[CHEMFP_ALIGN8_SMALL].method_p->intersect_popcount;
    }
  }

  /* Check for 4 byte alignment */

  if (ALIGNMENT(arena1, 4) == 0 &&
      ALIGNMENT(arena2, 4) == 0 &&
      storage_len1 % 4 == 0 &&
      storage_len2 % 4 == 0) {
    return alignments[CHEMFP_ALIGN4].method_p->intersect_popcount;
  }

  /* At least we're one byte aligned */
  return alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
}
