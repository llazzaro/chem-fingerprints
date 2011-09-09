#include "chemfp.h"
#include <stdio.h>

/* Bit operations related to byte and hex fingerprints

  A byte fingerprint is a length and a sequence of bytes where each byte
  stores 8 fingerprints bits, in the usual order. (That is, the byte 'A',
  which is the hex value 0x41, is the bit pattern "01000001".)


  A hex fingerprint is also stored as a length and a sequence of bytes
  but each byte encode 4 bits of the fingerprint as a hex character. The
  only valid byte values are 0-9, A-F and a-f. Other values will cause
  an error value to be returned. */

/***** Functions for hex fingerprints ******/

/* Map from ASCII value to bit count. Used with hex fingerprints.
   BIG is used in cumulative bitwise-or tests to check for non-hex input */

#define BIG 16
static int hex_to_value[256] = {
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  0,     1,   2,   3,   4,   5,   6,   7,   8,   9, BIG, BIG, BIG, BIG, BIG, BIG,

  /* Upper-case A-F */
  BIG,  10,  11,  12,  13,  14,  15, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,

  /* Lower-case a-f */
  BIG,  10,  11,  12,  13,  14,  15, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,

  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,

  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
};

/* Map from ASCII value to popcount. Used with hex fingerprints. */

static int hex_to_popcount[256] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   1,   1,   2,   1,   2,   2,   3,   1,   2,    0,   0,   0,   0,   0,   0,

    0,  2,   3,   2,   3,   3,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,  2,   3,   2,   3,   3,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,

    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,

    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
};

/* Map from an integer to its popcount. The maximum possible valid hex
   input is 'f'/'F', which is 15, but non-hex input will set bit 0x10, so
   I include the range 16-31 as well. */

int _popcount[32] = {
  0, 1, 1, 2, 1, 2, 2, 3,
  1, 2, 2, 3, 2, 3, 3, 4,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
};

/* Return 1 if the string contains only hex characters; 0 otherwise */
int chemfp_hex_isvalid(int len, const char *sfp) {
  int i, union_w=0;
  const unsigned char *fp = (unsigned char *) sfp;

  /* Out of range values set 0x10 so do cumulative bitwise-or and see if that
     bit is set. Optimize for the expected common case of validfingerprints. */
  for (i=0; i<len; i++) {
    union_w |= hex_to_value[fp[i]];
  }
  return (union_w < BIG) ? 1 : 0;
}

/* Return the population count of a hex fingerprint, otherwise return -1 */
int chemfp_hex_popcount(int len, const char *sfp) {
  int i, union_w=0, popcount=0;
  const unsigned char *fp = (const unsigned char *) sfp;

  for (i=0; i<len; i++) {
    /* Keep track of the cumulative popcount and the cumulative bitwise-or */
    popcount += hex_to_popcount[fp[i]];
    union_w |= hex_to_value[fp[i]];
  }
  if (union_w >= BIG) {
    return -1;  /* Then this was an invalid fingerprint (contained non-hex characters) */
  }
  return popcount;
}

/* Return the population count of the intersection of two hex fingerprints,
   otherwise return -1. */
int chemfp_hex_intersect_popcount(int len, const char *sfp1, const char *sfp2) {
  int i, union_w=0, intersect_popcount=0;
  int w1, w2;
  const unsigned char *fp1 = (const unsigned char *) sfp1;
  const unsigned char *fp2 = (const unsigned char *) sfp2;

  for (i=0; i<len; i++) {
    /* Get the popcount for each hex value. (Or 0 for non-hex values.) */
    w1 = hex_to_value[fp1[i]];
    w2 = hex_to_value[fp2[i]];
    /* Cumulative bitwise-or to check for non-hex values  */
    union_w = union_w | (w1|w2);
    intersect_popcount = intersect_popcount + _popcount[w1 & w2];
  }
  if (union_w >= BIG) {
    return -1;
  }
  return intersect_popcount;
}

/* Return the Tanitoto between two hex fingerprints, or -1.0 for invalid fingerprints
   If neither fingerprint has any set bits then return 1.0 */
/* I spent a lot of time trying out different ways to optimize this code.
   This is quite fast, but feel free to point out better ways! */
double chemfp_hex_tanimoto(int len, const char *sfp1, const char *sfp2) {
  int i=0, union_w=0;
  int union_popcount=0, intersect_popcount=0;
  int w1, w2;
  int w3, w4;
  int upper_bound = len - (len%2);
  const unsigned char *fp1 = (const unsigned char *) sfp1;
  const unsigned char *fp2 = (const unsigned char *) sfp2;

  /* Hex fingerprints really should be even-length since two hex characters
     are used for a single fingerprint byte and all chemfp fingerprints must
     be a multiple of 8 bits. I'll allow odd-lengths since I don't see how
     that's a bad idea and I can see how some people will be confused by
     expecting odd lengths to work. More specifically, I was confused because
     I used some odd lengths in my tests. ;) */

  /* I'll process two characters at a time. Loop-unrolling was about 4% faster. */
  for (; i<upper_bound; i+=2) {
    w1 = hex_to_value[fp1[i]];
    w2 = hex_to_value[fp2[i]];
    w3 = hex_to_value[fp1[i+1]];
    w4 = hex_to_value[fp2[i+1]];
    /* Check for illegal characters */
    union_w |= (w1|w2|w3|w4);
    /* The largest possible index is w1|w2 = (16 | 15) == 31 and */
    /* is only possible when the input is not a legal hex character. */
    union_popcount += _popcount[w1|w2]+_popcount[w3|w4];
    /* The largest possible index is w1&w2 = (16 & 16) == 16 */
    intersect_popcount += _popcount[w1&w2]+_popcount[w3&w4];
  }
  /* Handle the final byte for the case of odd fingerprint length */
  for (; i<len; i++) {
    w1 = hex_to_value[fp1[i]];
    w2 = hex_to_value[fp2[i]];
    /* Check for illegal characters */
    union_w |= (w1|w2);
    /* The largest possible index is (16 | 15) == 31 */
    /* (and only when the input isn't a legal hex character) */
    union_popcount += _popcount[w1|w2];
    /* The largest possible index is (16 & 16) == 16 */
    intersect_popcount += _popcount[w1&w2];
  }
  /* Check for illegal character */
  if (union_w >= BIG) {
    return -1.0;
  }
  /* Special case define that 0/0 = 0.0. It's hard to decide what to 
	 use here, for example, OpenEye uses 1.0. It seems that 0.0
     is the least surprising choice. */
  if (union_popcount == 0) {
    return 0.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;  /* +0.0 to coerce to double */
}

/* Return 1 if the query fingerprint is contained in the target, 0 if it isn't,
   or -1 for invalid fingerprints */
/* This code assumes that 1) most tests fail and 2) most fingerprints are valid */
int chemfp_hex_contains(int len, const char *squery_fp,
                        const char *starget_fp) {
  int i, query_w, target_w;
  int union_w=0;
  const unsigned char *query_fp = (const unsigned char *) squery_fp;
  const unsigned char *target_fp = (const unsigned char *) starget_fp;

  for (i=0; i<len; i++) {
    /* Subset test is easy; check if query & target == query
       I'll do word-by-word tests, where the word can also overflow to BIG
       Do the normal test against BIG to see if there was a non-hex input */
    query_w = hex_to_value[query_fp[i]];
    target_w = hex_to_value[target_fp[i]];
    union_w |= (query_w|target_w);
    if ((query_w & target_w) != query_w) {
      /* Not a subset, but first, check if there was a a non-hex input */
      if (union_w >= BIG) {
        return -1;
      }
      return 0;
    }
  }
  /* This was a subset, but there might have been a non-hex input */
  if (union_w >= BIG) {
    return -1;
  }
  return 1;
}

/****** byte fingerprints *******/

/* These algorithms are a lot simpler than working with hex fingeprints.
   There are a number of performance tweaks I could put in, especially
   if I know the inputs are word aligned, but I'll leave those for later. */

static int byte_popcounts[] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8  };


/* Return the population count of a byte fingerprint */
/* There are faster algorithms but this one is fast and simple, and it doesn't
   place any requirements on word alignment. */
int chemfp_byte_popcount(int len, const unsigned char *fp) {
  int i, popcount = 0;
  for (i=0; i<len; i++) {
    popcount += byte_popcounts[fp[i]];
  }
  return popcount;
}

/* Return the population count of the intersection of two byte fingerprints */
int chemfp_byte_intersect_popcount(int len, const unsigned char *fp1,
                                   const unsigned char *fp2) {
  int i, intersect_popcount = 0;
  for (i=0; i<len; i++) {
    intersect_popcount += byte_popcounts[fp1[i]&fp2[i]];
  }
  return intersect_popcount;
}

/* Return the Tanitoto between two byte fingerprints, or -1.0 for invalid fingerprints
   If neither fingerprint has any set bits then return 1.0 */
double chemfp_byte_tanimoto(int len, const unsigned char *fp1,
                            const unsigned char *fp2) {
  int i, union_popcount=0, intersect_popcount=0;
  /* Accumulate the total union and intersection popcounts */
  for (i=0; i<len; i++) {
    union_popcount += byte_popcounts[fp1[i] | fp2[i]];
    intersect_popcount += byte_popcounts[fp1[i] & fp2[i]];
  }
  /* Special case for when neither fingerprint has any bytes set */
  if (union_popcount == 0) {
    return 1.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;  /* +0.0 to coerce to double */
}

/* Return 1 if the query fingerprint is contained in the target, 0 if it isn't */
int chemfp_byte_contains(int len, const unsigned char *query_fp,
                         const unsigned char *target_fp) {
  int i;
  for (i=0; i<len; i++) {
    if ((query_fp[i] & target_fp[i]) != query_fp[i]) {
      return 0;
    }
  }
  return 1;
}


/* Return the Tanitoto between a byte fingerprint and a hex fingerprint */
/* The size is the number of bytes in the byte_fp */
double chemfp_byte_hex_tanimoto(int size,
				const unsigned char *byte_fp,
				const char *shex_fp) {
  const unsigned char *hex_fp = (unsigned char *) shex_fp;
  int union_w=0;
  int union_popcount=0, intersect_popcount=0;
  int w1, w2;
  int byte;
  unsigned char wc;

  /* I'll process two characters at a time. Loop-unrolling was about 4% faster. */
  while (size > 0) {
    w1 = hex_to_value[*hex_fp++];
    w2 = hex_to_value[*hex_fp++];
    /* Check for illegal characters */
    union_w |= (w1|w2);
    wc = (w1<<4) | w2;
    byte = *byte_fp++;
    union_popcount += byte_popcounts[byte | wc];
    intersect_popcount += byte_popcounts[byte & wc];
    size--;
  }
  if (union_w >= BIG) {
    return -1.0;
  }
  /* Special case define that 0/0 = 0.0. It's hard to decide what to 
	 use here, for example, OpenEye uses 1.0. It seems that 0.0
     is the least surprising choice. */
  if (union_popcount == 0) {
    return 0.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;  /* +0.0 to coerce to double */
}
