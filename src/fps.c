/* This file contains code which deals with the hex file format */

#include <stdio.h>
#include <string.h>

#include "heapq.h"
#include "chemfp.h"


/* Use <= to preserve order */
/* Python changed from <= to < . See comments in heapq.c */
/* Should I use the unique index as the secondary key? */
/* That would enforce ordering. */
static int fps_heap_lt(chemfp_heap *heap, int i, int j) {
  if (heap->scores[i] < heap->scores[j])
	return 1;
  if (heap->scores[i] > heap->scores[j])
	return 0;
  // break ties on a first-come basis
  return (heap->indicies[i] > heap->indicies[j]);
}

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


// newline must be known to exist
// This will stop at the newline.
// That's why there's no need to have the line length
static int parse_line(int hex_len, char *line,
					  char **id_start, int *id_len) {
  int fp_field_len, ws_len, tmp_id_len;
  char *s;

  fp_field_len = strspn(line, "0123456789abcdefABCDEF");
  if (fp_field_len == 0)
	return CHEMFP_MISSING_FINGERPRINT;
  if (fp_field_len % 2 != 0)
	return CHEMFP_BAD_FINGERPRINT;
  if (hex_len != -1 && hex_len != fp_field_len)
	return CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH;

  s = line+fp_field_len;
  // What if there's a NUL there? That's not whitespace, so it
  // will be part of the fingerprint. But that's not allowed.
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

  // This stops as the next whitespace. There must be a whitespace
  // because I checked that the block ends with a newline.
  // Also check for illegal whitespace.
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

static char *chemfp_to_next_line(char *s) {
  while (*s != '\n')
	s++;
  return s+1;
}


int chemfp_fps_line_validate(int hex_len, int line_len, char *line) {
  char *id_start;
  int id_end;
  // check for obviously bad parameters? That's not the C way.
  //if (len<0)
  //  return CHEMFP_BAD_ARG;
  // However, checking user input is essential.
  if (line_len==0 || line[line_len-1] != '\n')
	return CHEMFP_MISSING_NEWLINE;
  return parse_line(hex_len, line, &id_start, &id_end);
}

int chemfp_fps_tanimoto(int hex_len, char *hex_query,
						int target_block_len, char *target_block,
						double threshold,
						int *num_found_p,
						char **id_starts, int *id_lens,
						double *scores,
						int *lineno_p) {
  int lineno = (lineno_p ? *lineno_p : 1);
  char *line = target_block;
  char *end = target_block + target_block_len;
  int num_found = 0, err;
  char *s;
  double score;

  if (target_block_len==0 || target_block[target_block_len-1] != '\n') {
	*num_found_p = -1; // XXX right?
	return CHEMFP_MISSING_NEWLINE;
  }
  while (line < end) {
	err = parse_line(hex_len, line, id_starts, id_lens);
	if (err < 0)
	  goto finish;

	// Find the end of the line.
	// It might already be at the end of the line.
	s = chemfp_to_next_line(*id_starts + *id_lens);

	score = chemfp_hex_tanimoto(hex_len, hex_query, line);
	if (score >= threshold) {
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
	*num_found_p = -1; // XXX right?
	return CHEMFP_MISSING_NEWLINE;
  }
  while (line < end) {
	err = parse_line(hex_len, line, &id_start, &id_len);
	if (err < 0)
	  goto finish;
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

  int size = heap->size;
  int k = heap->k;
  int unique_idx = heap->unique_idx;

  // Test if we're growing the heap to k, or if the heap is already at k
  if (size == k)
	goto replace_in_heap;

  while (line < end) {
	err = parse_line(hex_len, line, &id_start, &id_len);
	if (err < 0) {
	  heap->size = size;
	  goto finish;
	}

	next_line = chemfp_to_next_line(id_start+id_len);
	score = chemfp_hex_tanimoto(hex_len, hex_query, line);
	if (score >= threshold) {
	  heap->indicies[size] = unique_idx++;
	  heap->scores[size] = score;
	  heap->id_starts[size] = id_start;
	  heap->id_lens[size] = id_len;
	  size++;

	  if (size == k) {
		heapq_heapify(k, (void *)heap,
					  (heapq_lt) fps_heap_lt, (heapq_swap) fps_heap_swap);
		// Update the minimum threshold
		heap->threshold = threshold;
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
	  heapq_siftup(heap->size, (void *) heap, 0,
				   (heapq_lt) fps_heap_lt, (heapq_swap) fps_heap_swap);
	  threshold = heap->scores[0];
	  /*
		// could jump to the end here
		// But then the lineno isn't right for the block
		// and I don't do validation.
	  if (threshold == 1.0) {
	  }
	  */
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
	heapq_heapify(heap->size, (void *)heap,
				  (heapq_lt) fps_heap_lt, (heapq_swap) fps_heap_swap);
  }
  heapq_heapsort(heap->size, (void *)heap,
				 (heapq_lt) fps_heap_lt, (heapq_swap) fps_heap_swap);
}
