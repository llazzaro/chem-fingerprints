#include "chemfp.h"
#include <stdio.h>

/* There are faster implementation than these, but these */
/* ones are easy to understand and pretty fast. */


#define BIG 16
static int hex_to_value[256] = {
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  0,     1,   2,   3,   4,   5,   6,   7,   8,   9, BIG, BIG, BIG, BIG, BIG, BIG,

  BIG,  10,  11,  12,  13,  14,  15, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
  BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG, BIG,
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

int _popcount[32] = {
  0, 1, 1, 2, 1, 2, 2, 3,
  1, 2, 2, 3, 2, 3, 3, 4,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
};


int chemfp_hex_isvalid(int len, const unsigned char *fp) {
  int i, union_w=0;
  for (i=0; i<len; i++) {
	union_w |= hex_to_value[fp[i]];
  }
  return (union_w < BIG);
}
  

int chemfp_hex_popcount(int len, const unsigned char *fp) {
  int i, union_w=0, popcount=0;
  for (i=0; i<len; i++) {
	union_w |= hex_to_value[fp[i]];
	popcount += hex_to_popcount[fp[i]];
  }
  if (union_w >= BIG) {
	return -1;
  }
  return popcount;
}

int chemfp_hex_intersect_popcount(int len, const unsigned char *fp1,
								  const unsigned char *fp2) {
  int i, union_w=0, intersect_popcount=0;
  int w1, w2;

  for (i=0; i<len; i++) {
	w1 = hex_to_value[fp1[i]];
	w2 = hex_to_value[fp2[i]];
	union_w = union_w | (w1|w2);
	intersect_popcount = intersect_popcount + _popcount[w1 & w2];
  }
  if (union_w >= BIG) {
	return -1;
  }
  return intersect_popcount;
}

double chemfp_hex_tanimoto(int len, const unsigned char *fp1,
						   const unsigned char *fp2) {
  int i=0, union_w=0;
  int union_popcount=0, intersect_popcount=0;
  int w1, w2;
  int w3, w4;
  int upper_bound = len - (len%2);
  // loop unrolling, take two at a time
  // This gives me 4% better timings for Tanimoto count
  for (; i<upper_bound; i+=2) {
	w1 = hex_to_value[fp1[i]];
	w2 = hex_to_value[fp2[i]];
	w3 = hex_to_value[fp1[i+1]];
	w4 = hex_to_value[fp2[i+1]];
	// Check for illegal characters
	union_w |= (w1|w2|w3|w4);
	// The largest possible index is (16 | 15) == 31
	// (and only when the input isn't a legal hex character)
	union_popcount += _popcount[w1|w2]+_popcount[w3|w4];
	// The largest possible index is (16 & 16) == 16
	intersect_popcount += _popcount[w1&w2]+_popcount[w3&w4];
  }
  // Should I allow odd-lengths for the hex fingerprints?
  for (; i<len; i++) {
	w1 = hex_to_value[fp1[i]];
	w2 = hex_to_value[fp2[i]];
	// Check for illegal characters
	union_w |= (w1|w2);
	// The largest possible index is (16 | 15) == 31
	// (and only when the input isn't a legal hex character)
	union_popcount += _popcount[w1|w2];
	// The largest possible index is (16 & 16) == 16
	intersect_popcount += _popcount[w1&w2];
  }
  if (union_w >= BIG) {
	return -1.0;
  }
  if (union_popcount == 0) {
	return 1.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;
}


int chemfp_hex_contains(int len, const unsigned char *query_fp,
						const unsigned char *target_fp) {
  int i, query_w, target_w;
  int union_w=0;
  for (i=0; i<len; i++) {
	query_w = hex_to_value[query_fp[i]];
	target_w = hex_to_value[target_fp[i]];
	union_w |= (query_w|target_w);
	if ((query_w & target_w) != query_w) {
	  if (union_w >= BIG) {
		return -1;
	  }
	  return 0;
	}
  }
  if (union_w >= BIG) {
	return -1;
  }
  return 1;
}

/****** byte *******/

static int byte_popcounts[] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8  };

int chemfp_byte_popcount(int len, const unsigned char *fp) {
  int i, popcount = 0;
  for (i=0; i<len; i++) {
	popcount += byte_popcounts[fp[i]];
  }
  return popcount;
}

int chemfp_byte_intersect_popcount(int len, const unsigned char *fp1,
								   const unsigned char *fp2) {
  int i, intersect_popcount = 0;
  for (i=0; i<len; i++) {
	intersect_popcount += byte_popcounts[fp1[i]&fp2[i]];
  }
  return intersect_popcount;
}


double chemfp_byte_tanimoto(int len, const unsigned char *fp1,
							const unsigned char *fp2) {
  int i, union_popcount=0, intersect_popcount=0;
  for (i=0; i<len; i++) {
	union_popcount += byte_popcounts[fp1[i] | fp2[i]];
	intersect_popcount += byte_popcounts[fp1[i] & fp2[i]];
  }
  if (union_popcount == 0) {
	return 1.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;
}

#if 0
// Some experimental code for testing a higher-performance popcount
typedef unsigned long long uint64;  //assume this gives 64-bits
typedef unsigned int uint32;  //assume this gives 32-bits
const uint64 m1  = 0x5555555555555555ULL; //binary: 0101...
const uint64 m2  = 0x3333333333333333ULL; //binary: 00110011..
const uint64 m4  = 0x0f0f0f0f0f0f0f0fULL; //binary:  4 zeros,  4 ones ...
const uint64 m8  = 0x00ff00ff00ff00ffULL; //binary:  8 zeros,  8 ones ...
const uint64 m16 = 0x0000ffff0000ffffULL; //binary: 16 zeros, 16 ones ...
const uint64 m32 = 0x00000000ffffffffULL; //binary: 32 zeros, 32 ones
const uint64 hff = 0xffffffffffffffffULL; //binary: all ones
const uint64 h01 = 0x0101010101010101ULL; //the sum of 256 to the power of 0,1,2,3...

// uintptr_t -- an int long enough for a pointer

// 1024 bits is 128 bytes is 16 uint64s
double chemfp_byte_tanimoto_128(const unsigned char *fp1,
								const unsigned char *fp2) {
  uint32 *ifp1 = (uint32 *) fp1;
  uint32 *ifp2 = (uint32 *) fp2;
  uint64 u, i;
  int union_popcount=0;
  int intersect_popcount=0;

  //  if ((fp1 & 0x7) != 4 || (fp2 & 0x7 != 4)) {
  //	printf("SKIP!\n");
  //	return chemfp_byte_tanimoto(128, fp1, fp2);
  //  }
  int n=16;
  do {
	u = (((uint64)(ifp1[0] | ifp2[0])) << 32) + (ifp1[1] | ifp2[1]);
	i = (((uint64)(ifp1[0] & ifp2[0])) << 32) + (ifp1[1] & ifp2[1]);
    u -= (u >> 1) & m1;             //put count of each 2 bits into those 2 bits
    i -= (i >> 1) & m1;
    u = (u & m2) + ((u >> 2) & m2); //put count of each 4 bits into those 4 bits 
    i = (i & m2) + ((i >> 2) & m2);
    u = (u + (u >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    i = (i + (i >> 4)) & m4;
    union_popcount += (u * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24)+...
    intersect_popcount += (i * h01)>>56;
	ifp1+=2;
	ifp2+=2;
  } while (--n);
  if (union_popcount == 0) {
	return 1.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;
}
#endif

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
