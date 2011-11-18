#ifndef CHEMFP_INTERNAL_H
#define CHEMFP_INTERNAL_H

#define ALIGNMENT(POINTER, BYTE_COUNT) \
  (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT))

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


#endif
