/* Functions for using the "fps" hex-based fingerprint file format */

#include <stdio.h>
#include <string.h>

#include "heapq.h"
#include "chemfp.h"

/* Internal function to parse a line from the fps block.
   The line MUST match /[0-9A-Fa-f]+\s\S+.*\n/
   It MUST be known to end with a newline. */
static int parse_line(
        int hex_len,  /* The expected length of the hex field, or -1 if unknown
                           (If it's known then it's used to validate.) */
        char *line,   /* The input line */
        char **id_start, int *id_len  /* After a successful return, these will contain
                                         the start position and length of the id field */
        ) {
  int fp_field_len, ws_len, tmp_id_len;
  char *s;

  /* Find the hex fingerprint and check that the length is appropriate */
  fp_field_len = strspn(line, "0123456789abcdefABCDEF");
  if (fp_field_len == 0)
    return CHEMFP_MISSING_FINGERPRINT;
  if (fp_field_len % 2 != 0)
    return CHEMFP_BAD_FINGERPRINT;
  if (hex_len != -1 && hex_len != fp_field_len)
    return CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH;

  s = line+fp_field_len;
  /* The only legal thing here is a space or a tab. */
  /* There might be some other character, including a NUL */
  /* XXX Why do I allow multiple whitespace ? Check the spec! */
  ws_len = strspn(s, " \t");  // \v? \f? \r?
  if (ws_len == 0) {
    switch (s[0]) {
    case '\n': return CHEMFP_MISSING_ID;
    case '\r': if (s[1] == '\n') return CHEMFP_MISSING_ID; /* else fallthrough */
    case '\v':
    case '\f': return CHEMFP_UNSUPPORTED_WHITESPACE;
    default: return CHEMFP_BAD_FINGERPRINT;
    }
  }
  s += ws_len;

  /* You must pass in a newline-terminated string to this function.
     Therefore, this function will finish while inside the string.
     Note that I'm also checking for illegal whitespace here. */
  tmp_id_len = strcspn(s, " \t\n\v\f\r");
  switch (s[tmp_id_len]) {
  case '\0': return CHEMFP_BAD_ID;
  case '\v':
  case '\f': return CHEMFP_UNSUPPORTED_WHITESPACE;
  case '\r': if (s[tmp_id_len+1] != '\n') return CHEMFP_UNSUPPORTED_WHITESPACE;
    break;
  }
  *id_start = s;
  *id_len = tmp_id_len;
  return CHEMFP_OK;
}

/* Go to the start of the next line. s may be at a newline already. */
static char *chemfp_to_next_line(char *s) {
  while (*s != '\n')
    s++;
  return s+1;
}


int chemfp_fps_line_validate(int hex_len, int line_len, char *line) {
  char *id_start;
  int id_end;
  if (line_len==0 || line[line_len-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;
  return parse_line(hex_len, line, &id_start, &id_end);
}

/* Compute Tanimoto scores for each line in the fps block and report all
   scores which are greater than or equal to the specified threshold. Callers
   must preallocate enough space in id_starts, id_lens, and scores for the
   results. */
int chemfp_fps_tanimoto(int hex_len, char *hex_query,
                        int target_block_len, char *target_block,
                        double threshold,
                        int *num_found_p,
                        char **id_starts, int *id_lens,
                        double *scores,
                        int *lineno_p) {
  /* The code is easier if I always track a line number even if not needed */
  int lineno = (lineno_p ? *lineno_p : 1);
  char *line = target_block;
  char *end = target_block + target_block_len;
  int num_found = 0, err;
  char *s;
  double score;

  if (target_block_len==0 || target_block[target_block_len-1] != '\n') {
    *num_found_p = 0;
    return CHEMFP_MISSING_NEWLINE;
  }
  while (line < end) {
    /* Parse a line, get the id start position and length, and verify hex_len */
    err = parse_line(hex_len, line, id_starts, id_lens);
    if (err < 0)
      goto finish;
    /* The character after the id might be a newline, or there might be other fields */
    s = chemfp_to_next_line(*id_starts + *id_lens);

    score = chemfp_hex_tanimoto(hex_len, hex_query, line);
    if (score >= threshold) { /* Record the match */
      *scores++ = score;
      id_starts++;
      id_lens++;
      num_found++;
    }
    line = s;
    lineno++;
  }
  err = CHEMFP_OK;
 finish:
  if (lineno_p)
    *lineno_p = lineno;
  *num_found_p = num_found;
  return err;
}

/* Return the number of fingerprints in the fps block which are greater
   than or equal to the specified threshold. */
int chemfp_fps_tanimoto_count(int hex_len, char *hex_query,
                              int target_block_len, char *target_block,
                              double threshold,
                              int *num_found_p, int *lineno_p) {
  int lineno = (lineno_p ? *lineno_p : 1);
  int num_found = *num_found_p;
  int id_len, err;
  char *line = target_block, *id_start, *s;
  char *end = target_block + target_block_len;
  double score;

  if (target_block_len == 0 || end[-1] != '\n') {
    *num_found_p = 0;
    return CHEMFP_MISSING_NEWLINE;
  }
  while (line < end) {
    /* Parse a line, get the id start position and length, and verify hex_len */
    err = parse_line(hex_len, line, &id_start, &id_len);
    if (err < 0)
      goto finish;
    /* The character after the id might be a newline, or there might be other fields */
    s = chemfp_to_next_line(id_start+id_len);

    score = chemfp_hex_tanimoto(hex_len, hex_query, line);
    if (score >= threshold)
      num_found++;
    line = s;
    lineno++;
  }
  err = CHEMFP_OK;
 finish:
  if (lineno_p)
    *lineno_p = lineno;
  *num_found_p = num_found;
  return err;
}

/****** Manage the best-of-N Tanimoto linear searches ********/

/* Compare two heap entries based on their score.
   Break ties based on the insertion index, with a preference to older entries. */
static int fps_heap_lt(chemfp_heap *heap, int i, int j) {
  if (heap->scores[i] < heap->scores[j])
    return 1;
  if (heap->scores[i] > heap->scores[j])
    return 0;
  // break ties on a first-come basis
  return (heap->indicies[i] > heap->indicies[j]);
}

/* Swap two entries in the heap */
static void fps_heap_swap(chemfp_heap *heap, int i, int j) {
  int idx = heap->indicies[i];
  double score = heap->scores[i];
  char *id_start = heap->id_starts[i];
  int id_len = heap->id_lens[i];

  heap->indicies[i] = heap->indicies[j];
  heap->scores[i] = heap->scores[j];
  heap->id_starts[i] = heap->id_starts[j];
  heap->id_lens[i] = heap->id_lens[j];

  heap->indicies[j] = idx;
  heap->scores[j] = score;
  heap->id_starts[j] = id_start;
  heap->id_lens[j] = id_len;
}


/* You have to set up memory space for everything. I just move some pointers around */
void chemfp_fps_heap_init(chemfp_heap *heap,
                          int k, double threshold,
                          int *indicies, double *scores,
                          char **id_starts, int *id_lens) {
  heap->size = 0;
  heap->k = k;
  heap->unique_idx = 0;
  heap->_reserved = 0;
  heap->threshold = threshold;
  heap->indicies = indicies;
  heap->scores = scores;
  heap->id_starts = id_starts;
  heap->id_lens = id_lens;
}

int chemfp_fps_heap_update_tanimoto(chemfp_heap *heap,
                                    int hex_len, char *hex_query,
                                    int target_block_len, char *target_block,
                                    int *lineno_p) {
  int lineno = (lineno_p ? *lineno_p : 1);
  char *line = target_block, *id_start;
  char *end = target_block + target_block_len;
  int err, id_len;
  char *next_line;
  double score;
  double threshold = heap->threshold;
  
  if (target_block_len == 0 || target_block[target_block_len-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;

  /* Data is pushed to this function, and processing state is stored
     in the heap. I have to figure out what I'm doing. */
  int size = heap->size;
  int k = heap->k;
  int unique_idx = heap->unique_idx;

  /* What I do is different if I've found k elements already. */
  if (size == k)
    goto replace_in_heap;

  /* Have not found k elements yet. The array values are still unsorted.  */
  while (line < end) {
    err = parse_line(hex_len, line, &id_start, &id_len);
    if (err < 0) {
      heap->size = size;  /* Save in case the error changed after already adding elements */
      goto finish;
    }

    next_line = chemfp_to_next_line(id_start+id_len);
    score = chemfp_hex_tanimoto(hex_len, hex_query, line);
    if (score >= threshold) {
      heap->indicies[size] = unique_idx++;  /* Append to the list */
      heap->scores[size] = score;
      heap->id_starts[size] = id_start;
      heap->id_lens[size] = id_len;
      size++;

      if (size == k) { /* Convert to a heap */
        chemfp_heapq_heapify(k, (void *)heap,
                 (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
        heap->threshold = threshold; /* "You must be at least -->this<-- tall to enter" */
        lineno++;
        heap->size = size;
        line = next_line;
        goto replace_in_heap;
      }
    }
    lineno++;
    line = next_line;
  }
  err = CHEMFP_OK;
  heap->size = size;
  goto finish;

 replace_in_heap:
  /* The heap already has k elements. New hits must replace old ones */
  while (line < end) {
    err = parse_line(hex_len, line, &id_start, &id_len);
    if (err < 0)
      goto finish;
    next_line = chemfp_to_next_line(id_start+id_len);
    score = chemfp_hex_tanimoto(hex_len, hex_query, line);
    if (score > threshold) {
      heap->indicies[0] = unique_idx++;
      heap->scores[0] = score;
      heap->id_starts[0] = id_start;
      heap->id_lens[0] = id_len;
      chemfp_heapq_siftup(heap->size, (void *) heap, 0,
                          (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
      threshold = heap->scores[0];
      if (threshold == 1.0) {
        /* No new element can be added because nothing can beat 1.0.
           However, I don't want to stop processing because this
           function also validates the fps block lines. I could go
           to another state, but that won't save any real time.
           So, I'll do nothing. */
      }
    }
    lineno++;
    line = next_line;
  }
  err = CHEMFP_OK;

 finish:
  heap->unique_idx = unique_idx;
  if (lineno_p)
    *lineno_p = lineno;
  return err;
}
void chemfp_fps_heap_finish_tanimoto(chemfp_heap *heap) {
  if (heap->size < heap->k) {
    // still adding to the array; need the first heapify
    chemfp_heapq_heapify(heap->size, (void *)heap,
                         (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
  }
  chemfp_heapq_heapsort(heap->size, (void *)heap,
                        (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
}
