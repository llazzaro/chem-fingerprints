#ifndef CPUID_H
#define CPUID_H

/**
 * @brief   Contains a portable cpuid implementation for x86 and
 *          x86-64 CPUs and POPCNT (SSE4.2) functions for molecular
 *          fingerprints.
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

#if defined(_MSC_VER) && (defined(_WIN32) || defined(_WIN64))
  #include <intrin.h>    /* __cpuid() */
  #include <nmmintrin.h> /* _mm_popcnt_u32(), _mm_popcnt_u64() */
#endif

#include <stdint.h>

/* %ecx bit flags */
#define bit_SSE3    (1 <<  0)
#define bit_SSSE3   (1 <<  9)
#define bit_SSE4_1  (1 << 19)
#define bit_SSE4_2  (1 << 20)
#define bit_POPCNT  (1 << 23)
#define bit_AVX     (1 << 28)

/* %edx bit flags */
#define bit_SSE     (1 << 25)
#define bit_SSE2    (1 << 26)

/**
 * Portable cpuid implementation for x86 and x86-64 CPUs
 * (supports PIC and non-PIC code).
 */
static void cpuid(unsigned int info,
                  unsigned int *eax,
                  unsigned int *ebx,
                  unsigned int *ecx,
                  unsigned int *edx)
{
#if defined(_MSC_VER) && (defined(_WIN32) || defined(_WIN64))
  int regs[4];
  __cpuid(regs, info);
  *eax = regs[0];
  *ebx = regs[1];
  *ecx = regs[2];
  *edx = regs[3];
#elif defined(__i386__) || defined(__i386)
  *eax = info;
  #if defined(__PIC__)
  __asm__ __volatile__ (
   "mov %%ebx, %%esi;" /* save %ebx PIC register */
   "cpuid;"
   "xchg %%ebx, %%esi;"
   : "+a" (*eax), 
     "=S" (*ebx),
     "=c" (*ecx),
     "=d" (*edx));
  #else
  __asm__ __volatile__ (
   "cpuid;"
   : "+a" (*eax), 
     "=b" (*ebx),
     "=c" (*ecx),
     "=d" (*edx));
  #endif
#elif defined(__x86_64__)
  *eax = info;
  __asm__ __volatile__ (
   "cpuid;"
   : "+a" (*eax), 
     "=b" (*ebx),
     "=c" (*ecx),
     "=d" (*edx));
#endif
}

/**
 * @return  An int value with the SSE and AVX bit flags set if the CPU
 *          supports the corresponding instruction sets.
 */
static int get_cpuid_flags()
{
  unsigned int info = 0x00000001;
  unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
  cpuid(info, &eax, &ebx, &ecx, &edx);
  return (edx & (bit_SSE    | 
                 bit_SSE2)) | 
         (ecx & (bit_SSE3   | 
                 bit_SSSE3  | 
                 bit_SSE4_1 | 
                 bit_SSE4_2 | 
                 bit_POPCNT | 
                 bit_AVX));
}

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
static int popcount_POPCNT(int size, const uint64_t *fp) {
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
static int intersect_popcount_POPCNT(int size, const uint64_t *fp1, const uint64_t *fp2) {
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

#endif /* CPUID_H */
