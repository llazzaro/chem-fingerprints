#include "popcount.h"

/* Quoting from Knuth, Fascicle 1,

  The first textbook on programming, "The Preparation of Programs for
  an Electronic Digital Computer" by Wilkes, Wheeler, and Gill,
  second edition (Reading, Mass.: Addison-Wesley, 1957), 155, 191-193,
  presented an interesting subroutine for sideways addition due to
  D. B. Gillies and J. C. P. Miller. Their method was devised for the
  35-bit numbers of the EDSAC, but it is readily converted to the
  following 64-bit procedure...

What follows is essentially this code, which is in Wikipedia
   http://en.wikipedia.org/wiki/Hamming_weight
as "popcount_3".

*/

int 
chemfp_popcount_gillies(int n, uint64_t *fp) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t h01 = UINT64_C(0x0101010101010101);

  int bit_count = 0, i;
  int size = (n+7) / 8;
  uint64_t x;

  for (i=0; i<size; i++) {
    x = fp[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}

int
chemfp_intersect_popcount_gillies(int n, uint64_t *fp1, uint64_t *fp2) {
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t h01 = UINT64_C(0x0101010101010101);

  int bit_count = 0, i;
  int size = (n+7) / 8;
  uint64_t x;

  for (i=0; i<size; i++) {
    x = fp1[i] & fp2[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bit_count += (int) ((x * h01) >> 56);
  }
  return bit_count;
}
