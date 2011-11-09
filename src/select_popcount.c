#include <stdio.h>
#include <stdint.h>

#include "chemfp.h"

/*#define REGULAR*/

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


static inline int popcount_lut8(int n, uint32_t *fp) {
  int cnt=0;
  unsigned int i;
  n /= 4;
  do {
    i = *fp;
    cnt += lut[i&255];
    cnt += lut[i>>8&255];
    cnt += lut[i>>16&255];
    cnt += lut[i>>24];  /* APD: removed the unneeded &255 */
    fp++;
  } while(--n);
  return cnt;
}

static inline int intersect_popcount_lut8(int n, uint32_t *fp1, uint32_t *fp2) {
  int cnt=0;
  unsigned int i;
  n /= 4;
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



/* Modified from http://stackoverflow.com/questions/1898153/how-to-determine-if-memory-is-aligned-testing-for-alignment-not-aligning */

#define ALIGNMENT(POINTER, BYTE_COUNT) \
  (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT))


/* chemfp stores fingerprints as Python strings */
/* (This may change in the future; memmap, perhaps?) */
/* The Python payload is 4 byte aligned but not 8 byte aligned. */

chemfp_popcount_f
chemfp_select_popcount(int num_bits,
		       int storage_len, const unsigned char *arena) {

  uintptr_t alignment = ALIGNMENT(arena, 4);
#ifdef REGULAR
    return chemfp_byte_popcount;
#endif
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

  if (alignment1 != alignment2) {
    /* No fast support yet for mixed alignments */
    return chemfp_byte_intersect_popcount;
  }

  if (alignment1 != 0) {
    /* I could in theory optimize this as well. */
    /* That will have to wait until there's a need for it. */
    return chemfp_byte_intersect_popcount;
  }

  if (storage_len1 % 4 != 0 ||
      storage_len2 % 4 != 0) {
    /* While fp[0] is 4 byte aligned, this means fp[1] is *NOT* */
    return chemfp_byte_intersect_popcount;
  }
  
  /* Yay! Four byte alignment! */

  if (num_bits < 256) {
    /* Probably not worthwhile */
    return chemfp_byte_intersect_popcount;
  }

  
  return (chemfp_intersect_popcount_f) intersect_popcount_lut8;

}
