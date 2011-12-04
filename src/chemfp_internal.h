#ifndef CHEMFP_INTERNAL_H
#define CHEMFP_INTERNAL_H

#define ALIGNMENT(POINTER, BYTE_COUNT) \
  (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT))


/* Macro to use for variable names which exist as a */
/* function parameter but otherwise aren't used */
/* This is to prevent compiler warnings on msvc /W4 */
#define UNUSED(x) (void)(x);


int _chemfp_get_option_report_popcount(void);
int _chemfp_set_option_report_popcount(int);

int _chemfp_get_option_report_intersect_popcount(void);
int _chemfp_set_option_report_intersect_popcount(int);

#endif
