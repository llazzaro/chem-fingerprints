#include "popcount.h"

#if defined(HAS_POPCOUNT_INTRINSIC)

int
_chemfp_popcount_intrinsic32(int n, uint32_t *fp) {
  int cnt=0;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    cnt += POPCNT32(*fp++);
  } while(--n);
  return cnt;
}

int
_chemfp_intersect_popcount_intrinsic32(int n, uint32_t *fp1, uint32_t *fp2) {
  int cnt=0;
  uint32_t i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+3) / 4;
  do {
    i = (*fp1 & *fp2);
    cnt += POPCNT32(i);
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}


int
_chemfp_popcount_intrinsic64(int n, uint64_t *fp) {
  int cnt=0;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+7) / 8;
  do {
    cnt += POPCNT64(*fp++);
  } while(--n);
  return cnt;
}

int
_chemfp_intersect_popcount_intrinsic64(int n, uint64_t *fp1, uint64_t *fp2) {
  int cnt=0;
  uint64_t i;

  /* Handle even cases where the fingerprint length is not a multiple of 4 */
  n = (n+7) / 8;
  do {
    i = (*fp1 & *fp2);
    cnt += POPCNT64(i);
    fp1++;
    fp2++;
  } while(--n);
  return cnt;
}

int
_chemfp_has_popcnt_instruction(void) {
  return 1;
}
#endif
