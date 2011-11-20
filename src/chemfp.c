#include "chemfp.h"

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
