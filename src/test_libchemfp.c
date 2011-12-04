#include <stdio.h>
#include <string.h>
#include "chemfp.h"

/* This is not a comprehensive test suite. That's done in Python code. */
/* This tests that the libchemfp public API is usable from C code */

#define CHECK(msg, expr, result)                                 \
  if ( (expr) != (result) ) {puts("FAIL: " msg); failed++;} \
  else {puts("PASS: " msg); passed++;}

int failed = 0;
int passed = 0;

int main() {
  char version_prefix[100];
  sprintf(version_prefix, "%d.%d", CHEMFP_MAJOR_VERSION, CHEMFP_MINOR_VERSION);

  puts("== info functions ==");
  CHECK("version", strncmp(chemfp_version(), version_prefix, strlen(version_prefix)), 0);
  CHECK("error code", strcmp(chemfp_strerror(CHEMFP_OK), "Ok"), 0);
  puts("== Hex checks ==");
  CHECK("empty string", chemfp_hex_isvalid(0, ""), 1);
  CHECK("2-bytes valid", chemfp_hex_isvalid(2, "abq"), 1);
  CHECK("all hex chars", chemfp_hex_isvalid(16, "0123456789abcdef"), 1);
  CHECK("invalid hex char", chemfp_hex_isvalid(16, "0123456789abcdeg"), 0);
  CHECK("popcount 0", chemfp_hex_popcount(4, "0000"), 0);
  CHECK("popcount 1", chemfp_hex_popcount(2, "01ff"), 1);
  CHECK("popcount bad", chemfp_hex_popcount(4, "01fg"), -1);
  CHECK("intersect popcount", chemfp_hex_intersect_popcount(6, "0F0123", "010b42"), 3);
  CHECK("Tanimoto empty", chemfp_hex_tanimoto(2, "00", "00"), 0.0);
  CHECK("Tanimoto", chemfp_hex_tanimoto(6, "123456", "012345"),
        (0.0+0+1+0+1+1)/(1+2+2+3+2+3));
  CHECK("Tanimoto fail", chemfp_hex_tanimoto(6, "12345 ", "012345"), -1.0);

  CHECK("contains", chemfp_hex_contains(2, "12", "3a"), 1);
  CHECK("does not contain", chemfp_hex_contains(2, "3a", "12"), 0);

  puts("== Binary checks ==");
  CHECK("empty popcount", chemfp_byte_popcount(0, "blah"), 0);
  CHECK("single byte popcount", chemfp_byte_popcount(1, "A"), 2);
  CHECK("multiple byte popcount", chemfp_byte_popcount(4, "ABCD"), 2+2+3+2);
  CHECK("intersect popcount", chemfp_byte_intersect_popcount(4, "ABCD", "BCDE"), 1+2+1+2);
  CHECK("Tanimoto null", chemfp_byte_tanimoto(0, "", ""), 0.0);
  CHECK("Tanimoto empty", chemfp_byte_tanimoto(1, "\0", "\0"), 0.0);
  CHECK("Tanimoto", chemfp_byte_tanimoto(2, "AB", "BC"), (1.0+2) / (3+3));
  CHECK("contains", chemfp_byte_contains(2, " *", "**"), 1);
  CHECK("does not contain", chemfp_byte_contains(2, "**", " *"), 0);
  
  puts("== FPS line ==");
  CHECK("valid with unknown hex len", chemfp_fps_line_validate(-1, 12, "abcdef\tspam\n"), 0);
  CHECK("valid with known hex len", chemfp_fps_line_validate(6, 12, "abcdef\tspam\n"), 0);
  CHECK("valid with bad hex len", chemfp_fps_line_validate(4, 12, "abcdef spam\n"),
        CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH);
  CHECK("valid with bad hex", chemfp_fps_line_validate(6, 12, "abcdeg spam\n"),
        CHEMFP_BAD_FINGERPRINT);

#if 0
  puts("== N-largest ==");
  int indices[2] = {0,0};
  double scores[2] = {0.0, 0.0};
  if (chemfp_nlargest_tanimoto_block(2,
                                     2, "A1",
                                     15, "her/hlvhSV#$(ZXVLzf*)4lkdf[]#@",
                                     0, -1,
                                     0.0,
                                     indices, scores) < 0) {
    puts("FAIL: chemfp_nlargest_tanimoto_block");
  } else {
    puts("PASS: chemfp_nlargest_tanimoto_block");
    //printf("[%d]=%f [%d]=%f\n", indices[0], scores[0], indices[1], scores[1]);
    CHECK("  indices[0]", indices[0], 10);
    CHECK("  indices[1]", indices[1], 13);
    CHECK("  scores[0]", scores[0], 0.375);
    CHECK("  scores[1]", scores[1], 0.363636363636363636);
  }
  
  scores[0] = scores[1] = -1;
  unsigned char *start_ids[2];
  int id_lens[2];
  int lineno;

  if (chemfp_hex_tanimoto_block(2,
                                4, "4131",
                                157,
"6865 ID0\n"
"722f ID1\n"
"686c ID2\n"
"7668 ID3\n"
"5356 ID4\n"
"2324 ID5\n"
"285a ID6\n"
"5856 ID7\n"
"4c7a ID8 blah\n"
"662a ID9\n"
"2934 ID10\n"
"6c6b ID11\n"
"6466 ID12\n"
"5b5d ID13 extra items\n"
"2340 ID14\n",
                                0.0,
                                scores, start_ids, id_lens, &lineno) < 0) {
    puts("FAIL: chemfp_hex_tanimoto_block");
  } else {
    puts("PASS: chemfp_hex_tanimoto_block");
    CHECK("  scores[0]", scores[0], 0.375);
    CHECK("  scores[1]", scores[1], 0.363636363636363636);
    CHECK("  id[0]", strncmp(start_ids[0], "ID10\n", id_lens[0]+1), 0);
    CHECK("  id[1]", strncmp(start_ids[1], "ID13 ", id_lens[1]+1), 0);
  }
  CHECK("chemfp_byte_intersect_popcount_count",
                // We know the BCDE returns an overlap of 6 bits
                // "A   " returns an overlap of 2
                chemfp_byte_intersect_popcount_count(4, "ABCD", 2, "XBCDEXA   X", 1, 5, 2),
                2);
  CHECK("chemfp_byte_intersect_popcount_count",
                chemfp_byte_intersect_popcount_count(4, "ABCD", 2, "XBCDEXA   X", 1, 5, 3),
                1);
#endif

  printf("Pass: %d   Fail: %d\n", passed, failed);
}
