#ifndef CPUID_H
#define CPUID_H

/**
 * @brief   Contains a portable cpuid implementation for x86 and
 *          x86-64 CPUs.
 * @author  Kim Walisch, <kim.walisch@gmail.com>
 * @version 1.1
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
  #include <intrin.h> /* __cpuid() */
#endif

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
 * @return  1 if the CPU supports the cpuid instruction else -1.
 */
static int cpuid(unsigned int info,
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
  return 1;
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
  return 1;
#elif defined(__x86_64__)
  *eax = info;
  __asm__ __volatile__ (
   "cpuid;"
   : "+a" (*eax), 
     "=b" (*ebx),
     "=c" (*ecx),
     "=d" (*edx));
  return 1;
#else
  /* compiler or CPU architecture do not support cpuid. */
  UNUSED(info);
  UNUSED(eax);
  UNUSED(ebx);
  UNUSED(ecx);
  UNUSED(edx);
  return -1;
#endif
}

/**
 * @return  An int value with the SSE and AVX bit flags set if the CPU
 *          supports the corresponding instruction sets.
 */
static int get_cpuid_flags(void)
{
  int flags = 0;
  unsigned int info = 0x00000001;
  unsigned int eax, ebx, ecx, edx;
  if (cpuid(info, &eax, &ebx, &ecx, &edx) != -1) {
    flags = (edx & (bit_SSE    | 
                    bit_SSE2)) | 
            (ecx & (bit_SSE3   | 
                    bit_SSSE3  | 
                    bit_SSE4_1 | 
                    bit_SSE4_2 | 
                    bit_POPCNT | 
                    bit_AVX));
  }
  return flags;
}

#endif /* CPUID_H */
