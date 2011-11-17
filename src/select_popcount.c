#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"

//#define REGULAR

/* Example of a faster algorithm assuming 4-byte aligned data */

static int lut[] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8  };


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
popcount_lauradoux(const uint64_t *fp, int size) {
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



typedef struct {
  int alignment;
  int min_size;
  int id;
  const char *name;
  chemfp_popcount_f popcount;
  chemfp_intersect_popcount_f intersect_popcount;
} method_type;

static method_type methods [] = {
  {1, 1, CHEMFP_METHOD_LUT8, "LUT8", 
   chemfp_byte_intersect_popcount, chemfp_byte_intersect_popcount},

  {2, 2, CHEMFP_METHOD_LUT16, "LUT16",  /* To implement */
   chemfp_byte_intersect_popcount, chemfp_byte_intersect_popcount},

  {2, 2, CHEMFP_METHOD_INTRINSIC16, "Intrinsic16",  /* To implement */
   chemfp_byte_intersect_popcount, chemfp_byte_intersect_popcount},

  {4, 4, CHEMFP_METHOD_INTRINSIC32, "Intrinsic32",  /* To implement */
   chemfp_byte_intersect_popcount, chemfp_byte_intersect_popcount},

  {8, 8, CHEMFP_METHOD_INTRINSIC64, "Intrinsic64",  /* To implement */
   chemfp_byte_intersect_popcount, chemfp_byte_intersect_popcount},

  {8, 96, CHEMFP_METHOD_LAURADOUX, "Lauradoux",
   popcount_lauradoux, intersect_popcount_lauradoux},
};

int
chemfp_num_methods(void) {
  return sizeof(methods) / sizeof(method_type);
}

const char *
chemfp_method_name(enum chemfp_popcount_methods method) {
  if (method < 0 || method >= chemfp_num_methods()) {
    return NULL;
  }
  return methods[method].name;
}


typedef struct {
  int alignment;
  int id;
  const char *name;
  method_type *method_p;
} alignment_type;

static alignment_type alignments[] = {
  {1, CHEMFP_ALIGN1, "1 byte aligned", &methods[CHEMFP_METHOD_LUT8]},
  {2, CHEMFP_ALIGN2, "2 byte aligned", &methods[CHEMFP_METHOD_LUT16]},
  {4, CHEMFP_ALIGN4, "4 byte aligned", &methods[CHEMFP_METHOD_LUT16]},
  {8, CHEMFP_ALIGN8_SMALL, "8 byte aligned (<96 bytes)", &methods[CHEMFP_METHOD_LUT16]},
  {8, CHEMFP_ALIGN8_LARGE, "8 byte aligned (>=96 bytes)", &methods[CHEMFP_METHOD_LAURADOUX]},
};

int
chemfp_num_alignments(void) {
  return sizeof(alignments) / sizeof(alignment_type);
}

const char *
chemfp_alignment_label(enum chemfp_alignments alignment) {
  if (alignment < 0 || alignment >= chemfp_num_alignments()) {
    return NULL;
  }
  return alignments[alignment].name;
}

int 
chemfp_get_alignment_method(enum chemfp_alignments alignment) {
  if (alignment < 0 || alignment >= chemfp_num_alignments()) {
    return -1;
  }
  return alignments[alignment].method_p->id;
}

int
chemfp_set_alignment_method(enum chemfp_alignments alignment,
			    enum chemfp_popcount_methods method) {
  if (alignment < 0 || alignment >= chemfp_num_alignments()) {
    return -1;
  }
  if (method < 0 || method >= chemfp_num_methods()) {
    return -1;
  }
  if (methods[method].alignment > alignments[alignment].alignment) {
    return -1;
  }
  alignments[alignment].method_p = methods+method;
  return 0;
}

  

/**************************************/

/* chemfp stores fingerprints as Python strings */
/* (This may change in the future; memmap, perhaps?) */
/* The Python payload is 4 byte aligned but not 8 byte aligned. */

chemfp_popcount_f
chemfp_select_popcount(int num_bits,
		       int storage_len, const unsigned char *arena) {

  uintptr_t alignment = ALIGNMENT(arena, 4);

  int num_bytes = (num_bits+7)/8;

#ifdef REGULAR
  return chemfp_byte_popcount;
#endif
  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }
  if (num_bytes < 4) {
    /* Really? That's very small. There's no reason to optimize this case */
    return chemfp_byte_popcount;
  }

  if (alignment != 0) {
    /* I could in theory optimize this as well. */
    return chemfp_byte_popcount;
  }

  if (storage_len % 4 != 0) {
    /* While fp[0] is 4 byte aligned, this means fp[1] is *NOT* */
    return chemfp_byte_popcount;
  }

  /* Yay! Four byte alignment! */
  return (chemfp_popcount_f) popcount_lut8;
}



chemfp_intersect_popcount_f
chemfp_select_intersect_popcount(int num_bits,
				 int storage_len1, const unsigned char *arena1,
				 int storage_len2, const unsigned char *arena2) {

  uintptr_t alignment1 = ALIGNMENT(arena1, 4);
  uintptr_t alignment2 = ALIGNMENT(arena2, 4);
  int storage_len = (storage_len1 < storage_len2) ? storage_len1 : storage_len2;

  int num_bytes = (num_bits+7)/8;

#ifdef REGULAR
  return chemfp_byte_intersect_popcount;
#endif

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }
  if (num_bytes < 8) {
    return alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
  }

  if (alignment1 != alignment2) {
    /* No fast support yet for mixed alignments */
    return alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
  }

  if (alignment1 != 0) {
    /* I could in theory optimize this case as well. */
    /* That will have to wait until there's a need for it. */
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

  /* Check for 2 byte alignment */

  if (ALIGNMENT(arena1, 2) == 0 &&
      ALIGNMENT(arena2, 2) == 0 &&
      storage_len1 % 2 == 0 &&
      storage_len2 % 2 == 0) {
    return alignments[CHEMFP_ALIGN2].method_p->intersect_popcount;
  }

  /* At least we're one byte aligned */
  return alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
}
