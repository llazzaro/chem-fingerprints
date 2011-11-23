/**
 * Enable the UINT64_C(c) macro from <stdint.h>.
 */
#if !defined(__STDC_CONSTANT_MACROS)
#  define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>

#include "popcount.h"

// popcount_lauradoux.c, v 1.1
// Created by Kim Walisch on 11/10/11.
//
// Changes:
//
// The main outer loop of the two functions now processes
// 12*8 = 96 bytes per iteration (previously 240 bytes).
// This makes the functions more efficient for small fingerprints
// e.g. 881 bits.

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
int
_chemfp_popcount_lauradoux(int byte_size, const uint64_t *fp) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
#if defined(ORIGINAL)
  const uint64_t h01 = UINT64_C(0x0101010101010101);
#endif
  
  uint64_t count1, count2, half1, half2, acc;
#if defined(ORIGINAL)
  uint64_t x;
#endif
  int i, j;    
  int size = (byte_size + 7) / 8;
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

#if !defined(ORIGINAL)
  /* Finish things up with the CHEMFP_ALIGN8_SMALL method */
  bit_count += _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->popcount(
			byte_size - limit*8, (unsigned char *) fp);
#else
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
#endif
  return bit_count;
}

int
_chemfp_intersect_popcount_lauradoux(int byte_size,
				     const uint64_t *fp1, const uint64_t *fp2) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
#if defined(ORIGINAL)
  const uint64_t h01 = UINT64_C(0x0101010101010101);
#endif
  
  uint64_t count1, count2, half1, half2, acc;
#if defined(ORIGINAL)
  uint64_t x;
#endif
  int i, j;
  /* Adjust the input size because the size isn't necessarily a multiple of 8 */
  /* even though the storage area will be a multiple of 8. */
  int size = (byte_size + 7) / 8;
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

#if !defined(ORIGINAL)
  /* Finish things up with the CHEMFP_ALIGN8_SMALL method */
  /* In my test case with 2048 bits the time went from 15.5 to 12.6 seconds */
  bit_count += _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->intersect_popcount(
			byte_size - limit*8,
			(unsigned char *) fp1, (unsigned char *) fp2);
    
#else
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
#endif
  return bit_count;
}
