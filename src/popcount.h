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
};

typedef int (*chemfp_method_check_f)(void);

typedef struct {
  int detected_index;
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


int _chemfp_popcount_lut8_4(int n, uint32_t *fp);
int _chemfp_intersect_popcount_lut8_4(int n, uint32_t *fp1, uint32_t *fp2);

int _chemfp_popcount_lut16_4(int n, uint32_t *fp);
int _chemfp_intersect_popcount_lut16_4(int n, uint32_t *fp1, uint32_t *fp2);

int _chemfp_popcount_lauradoux(int size, const uint64_t *fp);
int _chemfp_intersect_popcount_lauradoux(int size, const uint64_t *fp1, const uint64_t *fp2);


#if defined(_MSC_VER)
#  include <intrin.h>

#define HAS_POPCOUNT_INTRINSIC
#define POPCNT32(i) __popcnt(i)
#define POPCNT64(i) __popcnt64(i)

#else if defined(__GNUC__) || defined(__llvm__)

#define HAS_POPCOUNT_INTRINSIC
#define POPCNT32(i) __builtin_popcountl(i)
#define POPCNT64(i) __builtin_popcountll(i)

#endif

#if defined(HAS_POPCOUNT_INTRINSIC)
int _chemfp_popcount_intrinsic32(int n, uint32_t *fp);
int _chemfp_intersect_popcount_intrinsic32(int n, uint32_t *fp1, uint32_t *fp2);

int _chemfp_popcount_intrinsic64(int n, uint64_t *fp);
int _chemfp_intersect_popcount_intrinsic64(int n, uint64_t *fp1, uint64_t *fp2);

int _chemfp_has_popcnt_instruction(void);
#endif



#endif
