/**
 * @brief   Contains fast and portable bit population count functions
 *          for molecular fingerprints.
 * @author  Kim Walisch, <kim.walisch@gmail.com>
 * @version 1.1
 * @date    2011
 */

/** Enable the UINT64_C(c) macro from <stdint.h>. */
#if !defined(__STDC_CONSTANT_MACROS)
#  define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>

#include "popcount.h"

/**
 * Count the number of 1 bits (population count) in a fingerprint
 * using 64-bit tree merging. This implementation uses only 8
 * operations per 8 bytes on 64-bit CPUs.
 * The algorithm is due to Cédric Lauradoux, it is described and
 * benchmarked against other bit population count solutions (lookup
 * tables, bit-slicing) in his paper:
 * http://perso.citi.insa-lyon.fr/claurado/ham/overview.pdf
 * http://perso.citi.insa-lyon.fr/claurado/hamming.html
 */
int
chemfp_popcount_lauradoux(int byte_size, const uint64_t *fp) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
  const uint64_t h01 = UINT64_C(0x0101010101010101);
  
  uint64_t count1, count2, half1, half2, acc;
  uint64_t x;
  int i, j;
  int size = (byte_size + 7) / 8;
  int limit = size - size % 12;
  int bit_count = 0;

  /* The main outer loop processes 12*8 = 96 bytes per iteration
     (previously 240 bytes). This makes the popcount more efficient
     for small fingerprints e.g. 881 bits. */
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

  /* Count the bits of the remaining bytes (MAX 88) using 
     "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
     the code uses wikipedia's 64-bit popcount_3() implementation:
     http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation */

  /* Note: This is the "Gillies" algorithm, and timing tests show that
     it's more effective to put it here than to call the method tied
     to CHEMFP_ALIGN8_SMALL. You might think it best to use the
     fastest algorithm, but if you are using Lauradoux then you are on
     a machine which doesn't have POPCNT and where Lauradoux is faster
     than the lookup table. That's the same type of machine where
     Gillies is also faster than a lookup table.

     Calling it here instead of through the the function table saves
     time. It's 0.5% faster on my Mac (with gcc) and 5% faster on a
     Windows box (with msvc 10). */

  for (i = 0; i < size - limit; i++) {
    x = fp[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}

int
chemfp_intersect_popcount_lauradoux(int byte_size,
                                     const uint64_t *fp1, const uint64_t *fp2) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t m8  = UINT64_C(0x00FF00FF00FF00FF);
  const uint64_t m16 = UINT64_C(0x0000FFFF0000FFFF);
  const uint64_t h01 = UINT64_C(0x0101010101010101);
  uint64_t count1, count2, half1, half2, acc;
  uint64_t x;
  int i, j;
  int size = (byte_size + 7) / 8;
  int limit = size - size % 12;
  int bit_count = 0;

  /* 64-bit tree merging */
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

  /* intersect count the bits of the remaining bytes (MAX 88) using 
     "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
     the code uses wikipedia's 64-bit popcount_3() implementation:
     http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation */
  for (i = 0; i < size - limit; i++) {
    x = fp1[i] & fp2[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}
