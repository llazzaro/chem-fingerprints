#ifndef CHEMFP_INTERNAL_H
#define CHEMFP_INTERNAL_H

#define ALIGNMENT(POINTER, BYTE_COUNT) \
  (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT))

int _chemfp_get_option_report_popcount(void);
int _chemfp_set_option_report_popcount(int);

int _chemfp_get_option_report_intersect_popcount(void);
int _chemfp_set_option_report_intersect_popcount(int);

#endif
