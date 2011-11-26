#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"
#include "popcount.h"
#include "cpuid.h"

static unsigned long 
timeit(chemfp_popcount_f popcount, int size, int repeat);

/* These are the alignment categories which I support */

chemfp_alignment_type _chemfp_alignments[] = {
  {"align1", 1, 1, NULL},
  {"align4", 4, 4, NULL},
  {"align8-small", 8, 8, NULL},
  {"align8-large", 8, 96, NULL},

  /* This is a purely hack category. It's only used if set to "ssse3" */
  {"align-ssse3", 64, 64, NULL},
};

static int
has_popcnt_instruction(void) {
  return (get_cpuid_flags() & bit_POPCNT);
}

static int
has_ssse3(void) {
  return (get_cpuid_flags() & bit_SSSE3);
}

/* These are in the same order as an enum in popcount.h */
static chemfp_method_type compile_time_methods[] = {
  {0, CHEMFP_LUT8_1, "LUT8-1", 1, 1, NULL,
   _chemfp_popcount_lut8_1, _chemfp_intersect_popcount_lut8_1},

  {0, CHEMFP_LUT8_4, "LUT8-4", 4, 4, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut8_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut8_4},

  {0, CHEMFP_LUT16_4, "LUT16-4", 4, 4, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lut16_4,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lut16_4},

  {0, CHEMFP_LAURADOUX, "Lauradoux", 8, 96, NULL,
   (chemfp_popcount_f) _chemfp_popcount_lauradoux,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_lauradoux},

  {0, CHEMFP_POPCNT, "POPCNT", 8, 8,
   has_popcnt_instruction,
   (chemfp_popcount_f) _chemfp_popcount_popcnt,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_popcnt},

  {0, CHEMFP_GILLIES, "Gillies", 8, 8, NULL,
   (chemfp_popcount_f) _chemfp_popcount_gillies,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_gillies},

  {0, CHEMFP_SSSE3, "ssse3", 64, 64, has_ssse3,
   (chemfp_popcount_f) _chemfp_popcount_SSSE3,
   (chemfp_intersect_popcount_f) _chemfp_intersect_popcount_SSSE3},

};


/* These are the methods which are actually available at run-time */
/* This list is used for the public API */

static chemfp_method_type *
detected_methods[sizeof(compile_time_methods)/sizeof(chemfp_method_type)];

static int num_methods = 0;

static void
detect_methods(void) {
  int i, j=0;
  if (num_methods != 0) {
    return;
  }
  /* Go through all of the compile-time methods and see if it's available */
  for (i=0; i<sizeof(compile_time_methods)/sizeof(chemfp_method_type); i++) {
    if ((compile_time_methods[i].check == NULL) ||
        (compile_time_methods[i].check())) {

      /* Add it to the list of detected methods, and tell it its index position */
      compile_time_methods[i].detected_index = j;
      detected_methods[j++] = &compile_time_methods[i];
    }
  }
  num_methods = j;
}



int
chemfp_get_num_methods(void) {
  detect_methods();
  return num_methods;
}

const char *
chemfp_get_method_name(int method) {
  if (method < 0 || method >= chemfp_get_num_methods()) {
    return NULL;
  }
  return detected_methods[method]->name;
}


static void
set_default_alignment_methods(void) {
  int lut_method, best64_method, large_method, ssse3_method;
  unsigned long first_time, lut8_time, lut16_time, lut_time;
  unsigned long gillies_time, best64_time, lauradoux_time;
#if defined(GENERATE_SSSE3)
  unsigned long ssse3_time;
#endif

  /* Make sure we haven't already initialized the alignments */
  if (_chemfp_alignments[0].method_p != NULL) {
    return;
  }

  /* Figure out which methods are available for this hardware */
  detect_methods();

  /* This is the only possibility for 1-byte aligned */
  _chemfp_alignments[CHEMFP_ALIGN1].method_p = &compile_time_methods[CHEMFP_LUT8_1];

  /* Now do some timing measurements and figure out which method is
     likely the fastest for this hardware. It's a bit tricky; consider
     what happens if a timeslice boundary happens while doing a
     test. I mostly fix that by doing the timing twice and using
     the fastest time.

     I could require everyone call chemfp_select_fastest_method,
     but this should be good enough for almost everyone. */


  /* For 4-byte aligned we use a LUT. */
  /* TODO: implement a POPCNT instruction-based method for 4-byte aligned code */
  /* (You really should use an 8-byte aligned arena in this case, so not a priority */

  /* On older hardware the LUT16 can be slower than the LUT8 */
  first_time = timeit(compile_time_methods[CHEMFP_LUT8_4].popcount, 128, 200);
  lut8_time = timeit(compile_time_methods[CHEMFP_LUT8_4].popcount, 128, 200);
  if (first_time < lut8_time) {
    lut8_time = first_time;
  }

  first_time = timeit(compile_time_methods[CHEMFP_LUT16_4].popcount, 128, 200);
  lut16_time = timeit(compile_time_methods[CHEMFP_LUT16_4].popcount, 128, 200);
  if (first_time < lut16_time) {
    lut16_time = first_time;
  }

  /* Which one is faster? */
  if (lut8_time < lut16_time) {
    lut_method = CHEMFP_LUT8_4;
    lut_time = lut8_time;
  } else {
    lut_method = CHEMFP_LUT16_4;
    lut_time = lut16_time;
  }

  _chemfp_alignments[CHEMFP_ALIGN4].method_p = &compile_time_methods[lut_method];

  /* Let's see if the Gillies method is faster */
  first_time = timeit(compile_time_methods[CHEMFP_GILLIES].popcount, 128, 200);
  gillies_time = timeit(compile_time_methods[CHEMFP_GILLIES].popcount, 128, 200);
  if (first_time < gillies_time) {
    gillies_time = first_time;
  }

  /* For 8-byte aligned code we always want to use the POPCNT instruction if it exists */
  if (has_popcnt_instruction()) {
    _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p = 
      _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p = 
      _chemfp_alignments[CHEMFP_ALIGN_SSSE3].method_p = 
      &compile_time_methods[CHEMFP_POPCNT];
  } else {

    /* No POPCNT? Then either the LUT or Gillies for the small case, */
    /* and perhaps Lauradoux for the large case */
    if (lut_time < gillies_time) {
      best64_time = lut_time;
      best64_method = lut_method;
    } else {
      best64_time = gillies_time;
      best64_method = CHEMFP_GILLIES;
    }

    _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p = &compile_time_methods[best64_method];

    first_time = timeit(compile_time_methods[CHEMFP_LAURADOUX].popcount, 128, 200);
    lauradoux_time = timeit(compile_time_methods[CHEMFP_LAURADOUX].popcount, 128, 200);
    if (first_time < lauradoux_time) {
      lauradoux_time = first_time;
    }

    if (lauradoux_time < best64_time) {
      large_method = CHEMFP_LAURADOUX;
      best64_time = lauradoux_time;
    } else {
      large_method = best64_method;
    }
    _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p = &compile_time_methods[large_method];


    ssse3_method = large_method;

#if defined(GENERATE_SSSE3)
    if (has_ssse3()) {
      first_time = timeit(compile_time_methods[CHEMFP_SSSE3].popcount, 128, 200);
      ssse3_time = timeit(compile_time_methods[CHEMFP_SSSE3].popcount, 128, 200);
      if (first_time < ssse3_time) {
        ssse3_time = first_time;
      }
      if (ssse3_time < best64_time) {
        ssse3_method = CHEMFP_SSSE3;
      }
    }
#endif
    _chemfp_alignments[CHEMFP_ALIGN_SSSE3].method_p = &compile_time_methods[ssse3_method];

  }
}


int
chemfp_get_num_alignments(void) {
  set_default_alignment_methods();
  return sizeof(_chemfp_alignments) / sizeof(chemfp_alignment_type);
}

const char *
chemfp_get_alignment_name(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return NULL;
  }
  return _chemfp_alignments[alignment].name;
}

int 
chemfp_get_alignment_method(int alignment) {
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return CHEMFP_BAD_ARG;
  }
  return _chemfp_alignments[alignment].method_p->detected_index;
}

int
chemfp_set_alignment_method(int alignment, int method) {
  /* Make sure it's an available alignment and method */
  if (alignment < 0 || alignment >= chemfp_get_num_alignments()) {
    return CHEMFP_BAD_ARG;
  }
  if (method < 0 || method >= chemfp_get_num_methods()) {
    return CHEMFP_BAD_ARG;
  }
  /* Make sure the alignment and sizes are good enough */
  if (detected_methods[method]->alignment > _chemfp_alignments[alignment].alignment) {
    return CHEMFP_METHOD_MISMATCH;
  }
  if (detected_methods[method]->min_size > _chemfp_alignments[alignment].min_size) {
    return CHEMFP_METHOD_MISMATCH;
  }
  _chemfp_alignments[alignment].method_p = detected_methods[method];
  return CHEMFP_OK;
}

  

/**************************************/

/* chemfp stores fingerprints as Python strings */
/* (This may change in the future; memmap, perhaps?) */
/* The Python payload is 4 byte aligned but not 8 byte aligned. */

chemfp_popcount_f
chemfp_select_popcount(int num_bits,
                       int storage_len, const unsigned char *arena) {

  int num_bytes = (num_bits+7)/8;

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_default_alignment_methods();

  if (num_bytes <= 1) {
    /* Really? */
    return _chemfp_alignments[CHEMFP_ALIGN1].method_p->popcount;
  }
  if (ALIGNMENT(arena, 8) == 0 &&
      storage_len % 8 == 0) {
    if (num_bytes >= 96) {
      return _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p->popcount;
    } else {
      return _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->popcount;
    }
  }

  if (ALIGNMENT(arena, 4) &&
      storage_len % 4 == 0) {
    return _chemfp_alignments[CHEMFP_ALIGN4].method_p->popcount;
  }

  return _chemfp_alignments[CHEMFP_ALIGN1].method_p->popcount;
}



chemfp_intersect_popcount_f
chemfp_select_intersect_popcount(int num_bits,
                                 int storage_len1, const unsigned char *arena1,
                                 int storage_len2, const unsigned char *arena2) {

  int storage_len = (storage_len1 < storage_len2) ? storage_len1 : storage_len2;
  int num_bytes = (num_bits+7)/8;

  if (num_bytes > storage_len) {
    /* Give me bad input, I'll give you worse output */
    return NULL;
  }

  set_default_alignment_methods();
  
  if (num_bytes <= 1) {
    return _chemfp_alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
  }

  /* Check for 8 byte alignment */

  if (ALIGNMENT(arena1, 8) == 0 &&
      ALIGNMENT(arena2, 8) == 0 &&
      storage_len1 % 8 == 0 &&
      storage_len2 % 8 == 0) {

    /* We only use SSSE3 if this alignment is identical to "CHEMFP_SSSE3" */
    if (_chemfp_alignments[CHEMFP_ALIGN_SSSE3].method_p->id == CHEMFP_SSSE3) {

      /* I'll try, but only if I have 64 byte alignment */
      if (ALIGNMENT(arena1, 64) == 0 &&
          ALIGNMENT(arena2, 64) == 0 &&
          storage_len1 % 64 == 0 &&
          storage_len2 % 64 == 0) {
        return _chemfp_alignments[CHEMFP_ALIGN_SSSE3].method_p->intersect_popcount;
      }
    }

    if (num_bytes >= 96) {
      return _chemfp_alignments[CHEMFP_ALIGN8_LARGE].method_p->intersect_popcount;
    } else {
      return _chemfp_alignments[CHEMFP_ALIGN8_SMALL].method_p->intersect_popcount;
    }
  }

  /* Check for 4 byte alignment */

  if (ALIGNMENT(arena1, 4) == 0 &&
      ALIGNMENT(arena2, 4) == 0 &&
      storage_len1 % 4 == 0 &&
      storage_len2 % 4 == 0) {
    return _chemfp_alignments[CHEMFP_ALIGN4].method_p->intersect_popcount;
  }

  /* At least we're one byte aligned */
  return _chemfp_alignments[CHEMFP_ALIGN1].method_p->intersect_popcount;
}


/*********** Automatically select the fastest method ***********/

#if defined(_MSC_VER)
  #include <windows.h> /* QueryPerformanceCounter(LARGE_INTEGER*) */
#else
  #include <sys/time.h>
#endif

static long long
high_resolution_timer(void) {
#if defined(_MSC_VER)
  LARGE_INTEGER counter;
  if (!QueryPerformanceCounter(&counter) || 
      counter.QuadPart == 0) {
    fprintf(stderr, "Error: high resolution timer not available!\n");
    return 0;
  }
  return counter.QuadPart;
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  /* return usecs */
  return tv.tv_sec*1000000+tv.tv_usec;
#endif
}

/* Use uint64_t so it's 64-bit/8 byte aligned */
/* The contents are randomly generated. */
/* My first version was too small, and caused the LUT to appear
   faster even when the Gillies was better on real data. */

static uint64_t popcount_buffer[256] = {
  0x9b649615d1a50133ull,
  0xf3b8dada0e8b43deull,
  0x0197e207e4b9af2bull,
  0x68a2ecc4053b1305ull,
  0x93d933ac2f41e28full,
  0xb460859e01b6f925ull,
  0xc2c1a9eacc9e4999ull,
  0xdc5237f8200aec07ull,
  0x9e3bbe45d6e67641ull,
  0xa49bed7d060407d4ull,
  0xcca5f2913af53c5bull,
  0xfdd53575aab7c21aull,
  0x76b82d57bfa5c9ddull,
  0x0d2a87ba7f2439edull,
  0x9ec6a4ee2a6999d4ull,
  0xb9ae55f1f402ac97ull,
  0x08bbc6d1719a56bdull,
  0x969e5ef023c9ed23ull,
  0x6b7f08af661a9db6ull,
  0xad394da52bbbe18dull,
  0xdf9c3e28aae1c460ull,
  0xcf82e77d4f02f1efull,
  0x1fb88cdb648008ecull,
  0xc7a2ab7ecb8f84f5ull,
  0xbf8ef6833f18d407ull,
  0xb9c7eafdb4653fa2ull,
  0x90114b93b87a8a1dull,
  0x6e572c9e42e5061cull,
  0xb694ec549eeabc20ull,
  0xb362909621b9a2c8ull,
  0xcadab7b921d3cd0aull,
  0xd27f7aef7e2a0c6full,
  0xaf5d649ca1d2eefdull,
  0x6fc389a822e5769cull,
  0xdc849b5da5c5a101ull,
  0x3011e28954c71b98ull,
  0xecc6f2bb9b24b9d3ull,
  0x13d0974bbdbe16b5ull,
  0xb50625ca9f3348eeull,
  0x91a7462492f11cbbull,
  0x5fe0ca6928b55722ull,
  0xa5d89c3149133253ull,
  0x84645ec3c2cf4be6ull,
  0x22fd27c4b7981d9aull,
  0x3f9869fee13b43d9ull,
  0x0683208def61ce16ull,
  0x26f9fd185d31a581ull,
  0x837b1ded3af58f74ull,
  0x52e0246315b38ad7ull,
  0xbde27bb52d771b42ull,
  0x7fc2cb4428e33ee2ull,
  0xe3511d67a78fb94eull,
  0xeac2042d93f9d5f2ull,
  0xf987675f01562dd5ull,
  0x49f0250c27805c24ull,
  0xc331de3409aa714cull,
  0x9f3774691ac74fafull,
  0x167a091ad590c514ull,
  0xe4fbcf7d8f0f2008ull,
  0xfbc4b0cb233b04f6ull,
  0x960590126cce716aull,
  0x1dc1c707f6cc348dull,
  0x274b57e30bd6d6d3ull,
  0x67525306591d1746ull,
  0xf99163b382488844ull,
  0xe94f9bf47dfb0b16ull,
  0xcbb738584662cebbull,
  0x56ee87587103f7e5ull,
  0xcd8ff0352714830dull,
  0x624dd08f67e90c4bull,
  0xfff1f1b5b1f92417ull,
  0xcd9d4fb51b05e32bull,
  0x43c85c5a7a69cdc4ull,
  0xa27e72305a33c247ull,
  0xc40882a6813e08f1ull,
  0xad2b48e065ca1768ull,
  0x1ffa6c9616288e30ull,
  0xeb83e3323610ff2bull,
  0xb520d27b4f3a3273ull,
  0x15470f6c7346b910ull,
  0x3397c4c5b5e9bdc6ull,
  0x85f3179422591e54ull,
  0x86db696004af1781ull,
  0x22a9e51e871984beull,
  0x2de8e4cdd4652a1cull,
  0xe70ef696e037662aull,
  0xfc67e1f7083e10f0ull,
  0x945105f1c12fc00dull,
  0x4d169c35fc28ddebull,
  0x5522d55800e2b719ull,
  0x618040f560444bedull,
  0xff91b03867854f0bull,
  0x5ce1bfaf57be27d0ull,
  0x81752ce65cf5ba9eull,
  0x98e499fe7f0f365eull,
  0x5aa2bc888ad924bdull,
  0xae2de7838420c59bull,
  0x42cda0012ae00ff1ull,
  0x7620f99214e30e2full,
  0xa0be3f23a80f82ceull,
  0x420edefc42cedb09ull,
  0x80fe957c6a2817ffull,
  0x355174b6692ff140ull,
  0x47653e206352c78aull,
  0x808f7214b82d7c59ull,
  0x5dfcfe4144c253d4ull,
  0x4b918724a9084523ull,
  0x3e0608080fc35d1bull,
  0xf23cfdfd8c0b219eull,
  0x55bfd8597cdba8f5ull,
  0x269c25c3799d723cull,
  0x91e53b39bfdca5deull,
  0x02b04e9b8e52e823ull,
  0xc53fe276534e5317ull,
  0x18bd1dc656174acaull,
  0x0e5b4b3a13772eebull,
  0xa1943806fca56da6ull,
  0x04a5016c4c0be049ull,
  0x977ba238079e1e0cull,
  0x2df9dbcc4e036035ull,
  0x86adc435f1414d29ull,
  0x4402f529defe1868ull,
  0x03dbf44c63afc870ull,
  0xfbfe185f7297e08aull,
  0xe717fd0019ef65edull,
  0x7918c2b6e9275ba4ull,
  0x24f5ee4355f022b3ull,
  0xc0ba7a6be52fe0a4ull,
  0x685aabb6a61f00d8ull,
  0x3fa62a93e20e9372ull,
  0xc201d0ade1f15de7ull,
  0x28cb5915df8a4912ull,
  0x517843f1c3f9928full,
  0x4632606437902d9aull,
  0x82f853fb34d514b7ull,
  0x00464a29dcb32cbcull,
  0x84e1c0073eee811full,
  0x6eb2e2781ce72271ull,
  0xe3f40911bc8845e9ull,
  0xe6f2aacb1dd4d080ull,
  0xa87b1b15af61762full,
  0x810e66188c97dbeaull,
  0xdb919c39003db0d6ull,
  0x18452ccd19197178ull,
  0x5fe005b938986834ull,
  0xb179f1f3b113509full,
  0xea27088977c864c2ull,
  0x4e524739e812d35eull,
  0xf76f7a7d15cc08dbull,
  0xc0b9a7c0251f7f58ull,
  0x319d8eb2f9334c6dull,
  0x65db68328c2d2d4dull,
  0xc260bbf348039ee2ull,
  0xc692e00595613bffull,
  0x90fec8d4b374484dull,
  0x8ebd5b2ff1de52dfull,
  0xd3781952d5254631ull,
  0x84196d92f8852097ull,
  0xdc621b34a1763da6ull,
  0x0799e73b826efc26ull,
  0x098532b1f427cd10ull,
  0xfb2b0735121a374eull,
  0x9f8d3d10f5108176ull,
  0x57ee9d46db4529aaull,
  0x7c8db1c2e675c649ull,
  0x9d8e3388f3ef4382ull,
  0x639b5c10b29fc572ull,
  0x011f05e93ec9c4aeull,
  0xec28a9716fd3f5a1ull,
  0x837c0d205aefb577ull,
  0x0099fd93cadcb971ull,
  0xf29e78eae535df65ull,
  0x3c1ca48f330a6d1dull,
  0xb734f3c83f57de82ull,
  0x42f85b65c22dc638ull,
  0x0c50c85af7d3a601ull,
  0xea8ced5869fbe2fdull,
  0xb0cc396bfd86be6dull,
  0xb3ea7c3295866ef9ull,
  0x36cf28b306426badull,
  0x590de78ae5300681ull,
  0x41f4e16df296c0bcull,
  0xaad908beff6a93a9ull,
  0x909d243860e863d0ull,
  0x1d574b777f6e2725ull,
  0xacb7e3a9b94bb2b2ull,
  0x3b4d173db0b61bf6ull,
  0x4ccc5649c6c02c51ull,
  0x8d851d80b1a90638ull,
  0x6ca86fac5976ba0aull,
  0x09b49bdb4a58e177ull,
  0x7da8938aa92fe6b7ull,
  0x0f10d2d164ab5260ull,
  0x410822b41fff8a8eull,
  0x13d8dd389fe19217ull,
  0x0d6fcf685fdca839ull,
  0xae9965f4e51c9094ull,
  0x3cc74eabd4b3574aull,
  0x616a5f30b4a1e0a2ull,
  0x01c995c3cf9cde82ull,
  0x083e3df79ed6d08dull,
  0x50ca7def49e9be55ull,
  0x6827bee9c7b104adull,
  0xb09c88041e5a1480ull,
  0xd7d6b3f8a5fd79d2ull,
  0xe9a2a7562deb9cbbull,
  0xc6df657d5d037eaaull,
  0xa0513198d897cf1bull,
  0x941721727391ffbbull,
  0xdd65e39bef1199cbull,
  0x4e1129988fcc1a78ull,
  0x57d5274d4189e641ull,
  0xcd78a6383892a6c2ull,
  0x5380e97a1e588b36ull,
  0x4b153a04ed4f2d4cull,
  0x78c74fdda5d88d5full,
  0xa838c19ff3a05996ull,
  0x64a935bf0b55a732ull,
  0xa5727c5fee927c99ull,
  0x584c550d5f7af1d7ull,
  0x7b15564ed80dd58bull,
  0x42db540eda52029cull,
  0x78f64d45305d7f6full,
  0x8b549a03a9806568ull,
  0x6fa3c48b2b01ba66ull,
  0xc56ccbe0f05d1511ull,
  0x8adcd70ff4730081ull,
  0xf3f19cc845fd5b7aull,
  0x0936f92d55e55133ull,
  0xfda06bcd399ae365ull,
  0xde0c5052f3e158a4ull,
  0x58584d0c5e3b7dddull,
  0x3c3eb71846edfeb7ull,
  0xc1080e17c84266ffull,
  0xb25fd442e286d778ull,
  0x568605346b044740ull,
  0x54ffc2f936a972a2ull,
  0x366b795d073f062bull,
  0x206dadf277bbf8b4ull,
  0x916749a7cdf5e525ull,
  0x0afce12439536907ull,
  0x9fce50346e346701ull,
  0x562fe8ffc572a020ull,
  0xbac08aa15dc2f3f6ull,
  0x992aea3d03fb66a9ull,
  0x9e6a37740d285aafull,
  0x11dfb9a7b6b4424aull,
  0xe220772a626e2f9dull,
  0xae5c0a22b8ab8f2dull,
  0x11496ae8d4258860ull,
  0x6f3e74167f908fe6ull,
  0x622f3431103aef5dull,
  0x608584c6e190403dull,
  0xc8f7ec331fa3110cull,
  0x5ef7066f95c03fa1ull,
  0x48924db0f5d40254ull,
};

static unsigned long 
timeit(chemfp_popcount_f popcount, int size, int repeat) {
  long long t1, t2;
  unsigned char *start_buffer, *end_buffer, *fp;
  int i;
  if (size > sizeof(popcount_buffer)) {
    size = sizeof(popcount_buffer);
  }
  t1 = high_resolution_timer();
  start_buffer = (unsigned char *) popcount_buffer;
  end_buffer = start_buffer + sizeof(popcount_buffer);

  for (i=0; i<repeat; i++) {
    for (fp=start_buffer; fp<end_buffer; fp += size) {
      popcount(size, fp);
    }
  }
  t2 = high_resolution_timer();
  return (unsigned long)(t2-t1);
}



int
chemfp_select_fastest_method(int alignment, int repeat) {
  int method, best_method=-1, old_method;
  int probe_size;
  unsigned long dt;
  unsigned long first_time, best_time=0;
  chemfp_method_type *method_p=NULL;

  old_method = chemfp_get_alignment_method(alignment);
  if (old_method < 0) {
    return old_method;
  }

  /* NOTE: probe_size must evenly divide sizeof(popcount_buffer); */
  if (alignment == CHEMFP_ALIGN8_SMALL) {
    probe_size = 64; /* 512 bits; must be < 96 bytes  */
  } else {
    probe_size = 2048/8;
  }

  for (method=0; method<chemfp_get_num_methods(); method++) {

    /* See if I can use this method */
    if (chemfp_set_alignment_method(alignment, method) < 0) {
      continue;
    }
    method_p = _chemfp_alignments[alignment].method_p;

    /* Time the performance; do it twice in a timeslice happens in the middle  */
    first_time = timeit(method_p->popcount, probe_size, repeat);
    dt = timeit(method_p->popcount, probe_size, repeat);
    if (first_time < dt) {
      dt = first_time;
    }
    
    if (best_method == -1 || dt < best_time) {
      best_method = method;
      best_time = dt;
    }
  }
  if (best_method == -1) {
    /* Shouldn't happen, but I want to be on the safe side. */
    best_method = old_method;
  }
  chemfp_set_alignment_method(alignment, best_method);

  return best_method;
}
