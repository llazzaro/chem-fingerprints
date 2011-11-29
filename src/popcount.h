#if !defined(POPCOUNT_H)
#define POPCOUNT_H

#include <stdint.h>
#include "chemfp.h"
#include "chemfp_internal.h"


enum {
  CHEMFP_ALIGN1=0,
  CHEMFP_ALIGN4,
  CHEMFP_ALIGN8_SMALL,
  CHEMFP_ALIGN8_LARGE,
  CHEMFP_ALIGN_SSSE3
};

/* These are in the same order as compile_time_methods */
enum {
  CHEMFP_LUT8_1=0,
  CHEMFP_LUT8_4,
  CHEMFP_LUT16_4,
  CHEMFP_LAURADOUX,
  CHEMFP_POPCNT,
  CHEMFP_GILLIES,
  CHEMFP_SSSE3
};

typedef int (*chemfp_method_check_f)(void);

typedef struct {
  int detected_index;
  int id;
  const char *name;
  int alignment;
  int min_size;
  chemfp_method_check_f check;
  chemfp_popcount_f popcount;
  chemfp_intersect_popcount_f intersect_popcount;
} chemfp_method_type;

typedef struct {
  const char *name;
  int alignment;
  int min_size;
  chemfp_method_type *method_p;
} chemfp_alignment_type;


extern chemfp_alignment_type _chemfp_alignments[];

int _chemfp_popcount_lut8_1(int n, const unsigned char *fp);
int _chemfp_intersect_popcount_lut8_1(int n, const unsigned char *fp1, const unsigned char *fp2);

int _chemfp_popcount_lut8_4(int n, uint32_t *fp);
int _chemfp_intersect_popcount_lut8_4(int n, uint32_t *fp1, uint32_t *fp2);

int _chemfp_popcount_lut16_4(int n, uint32_t *fp);
int _chemfp_intersect_popcount_lut16_4(int n, uint32_t *fp1, uint32_t *fp2);

int _chemfp_popcount_gillies(int n, uint64_t *fp);
int _chemfp_intersect_popcount_gillies(int n, uint64_t *fp1, uint64_t *fp2);

int _chemfp_popcount_lauradoux(int size, const uint64_t *fp);
int _chemfp_intersect_popcount_lauradoux(int size, const uint64_t *fp1, const uint64_t *fp2);

int _chemfp_popcount_popcnt(int size, const uint64_t *fp);
int _chemfp_intersect_popcount_popcnt(int size, const uint64_t *fp1, const uint64_t *fp2);

int _chemfp_popcount_SSSE3(int, const unsigned*);
int _chemfp_intersect_popcount_SSSE3(int, const unsigned*, const unsigned*);

#endif
