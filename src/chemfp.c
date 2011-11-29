#include "chemfp.h"
#include "chemfp_internal.h"
#include <string.h>

const char *chemfp_version(void) {
  return CHEMFP_VERSION_STRING;
}

const char *chemfp_strerror(int err) {
  switch (err) {
  case CHEMFP_OK: return "Ok";
  case CHEMFP_BAD_ARG: return "Bad argument";
   
  case CHEMFP_NO_MEM: return "Cannot allocate memory";

  case CHEMFP_UNSUPPORTED_WHITESPACE: return "Unsupported whitespace";
  case CHEMFP_MISSING_FINGERPRINT: return "Missing fingerprint field";
  case CHEMFP_BAD_FINGERPRINT: return "Fingerprint field is in the wrong format";
  case CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH: return "Fingerprint is not the expected length";
  case CHEMFP_MISSING_ID: return "Missing id field";
  case CHEMFP_BAD_ID: return "Id field is in the wrong format";
  case CHEMFP_MISSING_NEWLINE: return "Line must end with a newline character";

  case CHEMFP_METHOD_MISMATCH: return "Mismatch between popcount method and alignment type";

  default: return "Unknown error";
  }
}

/* Code to handle getting and setting options */

typedef int (*get_option_f)(void);
typedef int (*set_option_f)(int);


typedef struct {
  const char *name;
  get_option_f getter;
  set_option_f setter;
} chemfp_option_type;

const chemfp_option_type chemfp_options[] = {
  {"report-popcount", _chemfp_get_option_report_popcount, _chemfp_set_option_report_popcount},
  {"report-intersect", _chemfp_get_option_report_intersect_popcount,
   _chemfp_set_option_report_intersect_popcount},
};

int
chemfp_get_num_options(void) {
  return sizeof(chemfp_options) / sizeof(chemfp_option_type);
}

const char *
chemfp_get_option_name(int i) {
  if (i < 0 || i >= chemfp_get_num_options()) {
    return NULL;
  }
  return chemfp_options[i].name;
}

int
chemfp_get_option(const char *option) {
  int i;
  for (i=0; i<chemfp_get_num_options(); i++) {
    if (!strcmp(chemfp_options[i].name, option)) {
      return chemfp_options[i].getter();
    }
  }
  return CHEMFP_BAD_ARG;
}

int
chemfp_set_option(const char *option, int value) {
  int i;
  for (i=0; i<chemfp_get_num_options(); i++) {
    if (!strcmp(chemfp_options[i].name, option)) {
      return chemfp_options[i].setter(value);
    }
  }
  return CHEMFP_BAD_ARG;
}
