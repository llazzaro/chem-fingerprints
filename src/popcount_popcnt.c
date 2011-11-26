/**
 * @brief   Contains portable popcount functions using the POPCNT
 *          (SSE4.2) instruction for molecular fingerprints.
 * @author  Kim Walisch, <kim.walisch@gmail.com>
 * @version 1.0
 * @date    2011
 *
 * The code within this file has been tested successfully with the
 * following compilers and operating systems:
 *
 * GNU GCC 4.4                    Linux i386 & x86-64 
 * LLVM clang 2.8,                Linux i386 & x86-64
 * Oracle Solaris Studio 12.2,    Linux i386 & x86-64
 * Intel C++ Composer XE 2011,    Linux i386 & x86-64, Windows 7 64-bit
 * GNU GCC 3.3,                   VMware Linux i386
 * GNU GCC 4.6,                   VMware Linux i386 & x86-64
 * Apple llvm-gcc-4.2,            Mac OS X 10.7
 * Apple clang version 3.0,       Mac OS X 10.7
 * Microsoft Visual Studio 2010,  Windows 7 64-bit
 * MinGW-w64 GCC 4.6,             Windows 7 64-bit
 */

#include "popcount.h"
#include <stdint.h>

#if defined(_MSC_VER) && (defined(_WIN32) || defined(_WIN64))
  #include <nmmintrin.h> /* _mm_popcnt_u32(), _mm_popcnt_u64() */
#endif

/** Convenience functions for the POPCNT instruction. */

#if defined(_MSC_VER) && defined(_WIN64)
static uint64_t POPCNT64(uint64_t x) {
  return _mm_popcnt_u64(x);
}
#elif defined(_MSC_VER) && defined(_WIN32)
static uint32_t POPCNT32(uint32_t x) {
  return _mm_popcnt_u32(x);
}
#elif defined(__x86_64__)
static uint64_t POPCNT64(uint64_t x) {
/* GNU GCC >= 4.2 supports the POPCNT instruction */
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
  __asm__ ("popcnt %1, %0" : "=r" (x) : "0" (x));
#endif
  return x;
}
#elif defined(__i386__) || defined(__i386)
static uint32_t POPCNT32(uint32_t x) {
/* GNU GCC >= 4.2 supports the POPCNT instruction */
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
  __asm__ ("popcnt %1, %0" : "=r" (x) : "0" (x));
#endif
  return x;
}
#endif

/**
 * Count the number of bits set in a fingerprint using the
 * the POPCNT (SSE4.2) instruction.
 * @warning  Use (get_cpuid_flags() & bit_POPCNT) to test if
 *           the CPU supports the POPCNT instruction.
 */
int _chemfp_popcount_popcnt(int size, const uint64_t *fp) {
  int bit_count = 0;
  int i;
#if defined(_WIN64) || defined(__x86_64__)
  size = (size + 7) / 8;
  for (i = 0; i < size; i++)
    bit_count += (int) POPCNT64(fp[i]);
#elif defined(_WIN32) || defined(__i386__) || defined(__i386)
  const uint32_t* fp_32 = (const uint32_t*) fp;
  size = (size + 3) / 4;
  for (i = 0; i < size; i++)
    bit_count += (int) POPCNT32(fp_32[i]);
#endif
  return bit_count;
}

/**
 * Count the number of bits set within the intersection of two
 * fingerprints using the POPCNT (SSE4.2) instruction.
 * @warning  Use (get_cpuid_flags() & bit_POPCNT) to test if
 *           the CPU supports the POPCNT instruction.
 */
int _chemfp_intersect_popcount_popcnt(int size, const uint64_t *fp1, const uint64_t *fp2) {
  int bit_count = 0;
  int i;
#if defined(_WIN64) || defined(__x86_64__)
  size = (size + 7) / 8;
  for (i = 0; i < size; i++)
    bit_count += (int) POPCNT64(fp1[i] & fp2[i]);
#elif defined(_WIN32) || defined(__i386__) || defined(__i386)
  const uint32_t* fp1_32 = (const uint32_t*) fp1;
  const uint32_t* fp2_32 = (const uint32_t*) fp2;
  size = (size + 3) / 4;
  for (i = 0; i < size; i++)
    bit_count += (int) POPCNT32(fp1_32[i] & fp2_32[i]);
#endif
  return bit_count;
}
