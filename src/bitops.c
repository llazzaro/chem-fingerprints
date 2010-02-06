#include "chemfp.h"


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
	union_w |= (w1|w2);
	intersect_popcount += _popcount[w1 & w2];
  }
  if (union_w >= BIG) {
	return -1;
  }
  return intersect_popcount;
}

double chemfp_hex_tanimoto(int len, const unsigned char *fp1,
						   const unsigned char *fp2) {
  int i, union_w=0;
  int union_popcount=0, intersect_popcount=0;
  int w1, w2;
  for (i=0; i<len; i++) {
	w1 = hex_to_value[fp1[i]];
	w2 = hex_to_value[fp2[i]];
	// Check for illegal characters
	union_w |= (w1|w2);
	// The largest possible index is (16 | 15) == 31
	// (and only when the input isn't a legal hex character)
	union_popcount += _popcount[w1|w2];
	// The largest possible index is (16 & 16) == 16
	intersect_popcount += _popcount[w1 & w2];
  }
  if (union_w >= BIG) {
	return -1.0;
  }
  if (union_popcount == 0) {
	return 0.0;
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
	return 0.0;
  }
  return (intersect_popcount + 0.0) / union_popcount;
}

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
