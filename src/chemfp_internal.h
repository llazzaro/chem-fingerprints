#ifndef CHEMFP_INTERNAL_H
#define CHEMFP_INTERNAL_H

#define ALIGNMENT(POINTER, BYTE_COUNT) \
  (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT))


#endif
