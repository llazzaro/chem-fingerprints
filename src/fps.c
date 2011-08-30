/* Functions for using the "fps" hex-based fingerprint file format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "heapq.h"
#include "chemfp.h"

enum {ADD_TO_HEAP, REPLACE_IN_HEAP, MAXED_OUT_HEAP};

/* Internal function to find the id field in an FPS line */
/* (Which means the fingerprint field is from line to *id_start-1 ) */
/* The line MUST match /^[0-9A-Fa-f]+\t\S+/ */
/* REQUIRED: the line MUST end with a newline (this is not checked) */
int chemfp_fps_find_id(
        int hex_size,  /* The expected length of the hex field, or -1 if unknown
                           (If it's known then it's used to validate.) */
        const char *line,   /* The input line */
        const char **id_start,  /* After a successful return, these will contain */
        const char **id_end     /* the start and end+1 position of the id field */
        ) {
  int fp_field_len, id_len;
  const char *s;

  /* Find the hex fingerprint and check that the length is appropriate */
  fp_field_len = strspn(line, "0123456789abcdefABCDEF");
  if (fp_field_len == 0)
    return CHEMFP_MISSING_FINGERPRINT;
  if (fp_field_len % 2 != 0)
    return CHEMFP_BAD_FINGERPRINT;
  if (hex_size != -1 && hex_size != fp_field_len)
    return CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH;

  s = line+fp_field_len;
  /* The only legal thing here is a tab. */
  /* Check if it's some other character, including a NUL */
  switch (s[0]) {
  case '\t': break;  /* The only legal option. Everything else improves the error code */
  case '\n': return CHEMFP_MISSING_ID;
  case '\r': if (s[1] == '\n') return CHEMFP_MISSING_ID; /* else fallthrough */
  case ' ':
  case '\v':
  case '\f': return CHEMFP_UNSUPPORTED_WHITESPACE;
  default: return CHEMFP_BAD_FINGERPRINT;
  }
  s++;

  /* You must pass in a newline-terminated string to this function.
     Therefore, this function will finish while inside the string.
     Note that I'm also checking for illegal whitespace here. */
  id_len = strcspn(s, " \t\n\v\f\r");
  switch (s[id_len]) {
  case '\0': return CHEMFP_BAD_ID;
  case ' ':
  case '\v':
  case '\f': return CHEMFP_UNSUPPORTED_WHITESPACE;
  case '\r': if (s[id_len+1] != '\n') return CHEMFP_UNSUPPORTED_WHITESPACE;
    break;
  }
  *id_start = s;
  *id_end = s+id_len;
  return CHEMFP_OK;
}

/* Go to the start of the next line. s may be at a newline already. */
static const char *chemfp_to_next_line(const char *s) {
  while (*s != '\n')
    s++;
  return s+1;
}


int chemfp_fps_line_validate(int hex_size, int line_size, const char *line_start) {
  const char *id_start, *id_end;
  if (line_size == 0 || line_start[line_size-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;
  return chemfp_fps_find_id(hex_size, line_start, &id_start, &id_end);
}

#if 0
/* Compute Tanimoto scores for each line in the fps block and report all
   scores which are greater than or equal to the specified threshold. Callers
   must preallocate enough space in id_starts, id_lens, and scores for the
   results. */
int chemfp_fps_tanimoto(int hex_size, char *hex_query,
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
    /* Parse a line, get the id start position and length, and verify hex_size */
    err = chemfp_fps_find_id(hex_size, line, id_starts, id_lens);
    if (err < 0)
      goto finish;
    /* The character after the id might be a newline, or there might be other fields */
    s = chemfp_to_next_line(*id_starts + *id_lens);

    score = chemfp_hex_tanimoto(hex_size, hex_query, line);
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
#endif

/* Return the number of fingerprints in the fps block which are greater
   than or equal to the specified threshold. */
int chemfp_fps_tanimoto_count(
	int num_bits,
	int query_storage_size,
	const unsigned char *query_arena, int query_start, int query_end,
	const char *target_block, int target_block_end,
        double threshold,
	int *counts, int *num_lines_processed) {
  const unsigned char *query_fp;
  const char *line, *next_line, *end;
  int fp_size = (num_bits+7)/8;
  int num_lines = 0, query_index;
  const char *id_start, *id_end;
  int err;
  double score;

  end = target_block + target_block_end;
  if (target_block_end == 0 || end[-1] != '\n') {
    err = CHEMFP_MISSING_NEWLINE;
    goto finish;
  }
  line = target_block;
  while (line < end) {
    /* Parse a line, get the id start position and length, and verify hex_size */
    err = chemfp_fps_find_id(fp_size*2, line, &id_start, &id_end);
    if (err < 0)
      goto finish;
    /* The character after the id might be a newline, or there might be other fields */
    next_line = chemfp_to_next_line(id_end);

    query_fp = query_arena + query_start * query_storage_size;
    for (query_index=query_start; query_index<query_end;
	 query_index++, query_fp += fp_size) {
      score = chemfp_byte_hex_tanimoto(fp_size, query_fp, line);
      if (score >= threshold)
	counts[query_index]++;
    }
    num_lines++;
    line = next_line;
  }
  err = CHEMFP_OK;
 finish:
  *num_lines_processed = num_lines;
  return err;
}

/****** Linear Tanimoto search with threshold and unlimited number of hits ********/

int chemfp_fps_threshold_tanimoto_search(
	int num_bits,
	int query_storage_size,
	const unsigned char *query_arena, int query_start, int query_end,
	
	const char *target_block, int target_block_end,
        double threshold,
	int num_cells, chemfp_tanimoto_cell *cells,
	const char ** stopped_at, int *num_lines_processed, int *num_cells_processed) {
  const char *line = target_block;
  const char *next_line;
  const char *end = target_block+target_block_end;
  const char *id_start, *id_end;
  const unsigned char *query_fp;
  chemfp_tanimoto_cell *current_cell;
  double score;
  int query_index;
  int num_lines = 0, num_queries;
  int err, retval;
  int fp_size = (num_bits+7)/8;

  current_cell = cells;
  if (query_start >= query_end) {
    retval = CHEMFP_OK;
    goto finish;
  }
  num_queries = query_end - query_start;
  if (end[-1] != '\n') {
    // There's no guarantee that the missing newline is on "stopped_at"
    // In the Python API there's no way to trigger this through normal code.
    retval = CHEMFP_MISSING_NEWLINE;
    goto finish;
  }

  while (line < end) {
    if (num_cells < num_queries) {
      goto success;
    }
    err = chemfp_fps_find_id(2*fp_size, line, &id_start, &id_end);
    if (err < 0) {
      retval = err;
      goto finish;
    }
    next_line = chemfp_to_next_line(id_end);

    query_fp = query_arena + query_start * query_storage_size;
    for (query_index=query_start; query_index<query_end;
	 query_index++, query_fp += fp_size) {
      score = chemfp_byte_hex_tanimoto(fp_size, query_fp, line);
      if (score >= threshold) {
	current_cell->score = score;
	current_cell->query_index = query_index;
	current_cell->id_start = id_start - target_block;
	current_cell->id_end = id_end - target_block;
	current_cell++;
	num_cells--;
      }
    }
    line = next_line;
    num_lines++;
  }
 success:
  retval = CHEMFP_OK;
  
 finish:
  *stopped_at = line;
  *num_lines_processed = num_lines;
  *num_cells_processed = current_cell - cells;
  return retval;
}

/****** Manage the best-of-N Tanimoto linear searches ********/

/* Compare two heap entries based on their score.
   Break ties based on the insertion index, with a preference to older entries. */
static int fps_heap_lt(chemfp_fps_heap *heap, int i, int j) {
  if (heap->scores[i] < heap->scores[j])
    return 1;
  if (heap->scores[i] > heap->scores[j])
    return 0;
  // break ties on a first-come basis
  return (heap->indicies[i] > heap->indicies[j]);
}

/* Swap two entries in the heap */
static void fps_heap_swap(chemfp_fps_heap *heap, int i, int j) {
  int idx = heap->indicies[i];
  double score = heap->scores[i];
  char *id = heap->ids[i];

  heap->indicies[i] = heap->indicies[j];
  heap->scores[i] = heap->scores[j];
  heap->ids[i] = heap->ids[j];

  heap->indicies[j] = idx;
  heap->scores[j] = score;
  heap->ids[j] = id;
}

#if 0
/* You have to set up memory space for everything. I just move some pointers around */
void chemfp_fps_heap_init(chemfp_fps_heap *heap,
                          int k, double threshold,
                          int *indicies, double *scores,
                          int *id_starts, int *id_ends) {
  heap->size = 0;
  heap->k = k;
  heap->unique_idx = 0;
  heap->_reserved = 0;
  heap->threshold = threshold;
  heap->indicies = indicies;
  heap->scores = scores;
  heap->id_starts = id_starts;
  heap->id_ends = id_ends;
}

int chemfp_fps_heap_update_tanimoto(chemfp_fsp_heap *heap,
                                    int hex_size, char *hex_query,
                                    int target_block_len, char *target_block,
                                    int *lineno_p) {
  int lineno = (lineno_p ? *lineno_p : 1);
  const char *line = target_block, *id_start, *id_end, *next_line;
  const char *end = target_block + target_block_len;
  int err;
  double score;
  double threshold = heap->threshold;
  int size, k, unique_idx;
  
  if (target_block_len == 0 || target_block[target_block_len-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;

  /* Data is pushed to this function, and processing state is stored
     in the heap. I have to figure out what I'm doing. */
  size = heap->size;
  k = heap->k;
  unique_idx = heap->unique_idx;

  /* What I do is different if I've found k elements already. */
  if (size == k)
    goto replace_in_heap;

  /* Have not found k elements yet. The array values are still unsorted.  */
  while (line < end) {
    err = chemfp_fps_find_id(hex_size, line, &id_start, &id_end);
    if (err < 0) {
      heap->size = size;  /* Save in case the error changed after already adding elements */
      goto finish;
    }

    next_line = chemfp_to_next_line(id_end);
    score = chemfp_hex_tanimoto(hex_size, hex_query, line);
    if (score >= threshold) {
      heap->indicies[size] = unique_idx++;  /* Append to the list */
      heap->scores[size] = score;
      heap->id_starts[size] = id_start-target_block;
      heap->id_ends[size] = id_end-target_block;
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
    err = chemfp_fps_find_id(hex_size, line, &id_start, &id_end);
    if (err < 0)
      goto finish;
    next_line = chemfp_to_next_line(id_end);
    score = chemfp_hex_tanimoto(hex_size, hex_query, line);
    if (score > threshold) {
      heap->indicies[0] = unique_idx++;
      heap->scores[0] = score;
      heap->id_starts[0] = id_start-target_block;
      heap->id_ends[0] = id_end-target_block;
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
#endif
/***************** new code */

int chemfp_fps_knearest_search_init(
        chemfp_fps_knearest_search *knearest_search,
	int num_bits, int query_storage_size,
	const unsigned char *query_arena, int query_start, int query_end,
        int k, double threshold) {

  chemfp_fps_heap *heaps = NULL;
  int *all_indicies = NULL;
  char **all_ids = NULL;
  double *all_scores = NULL;
  int i, num_queries;

  if (query_start >= query_end) {
    num_queries = 0;
    goto skip_malloc;
  } else {
    num_queries = query_end - query_start;
  }
  heaps = (chemfp_fps_heap *) calloc(sizeof(chemfp_fps_heap), num_queries);
  if (!heaps) {
    goto malloc_failure;
  }
  all_indicies = (int *) calloc(sizeof(int), k*num_queries);
  if (!all_indicies) {
    goto malloc_failure;
  }
  all_ids = (char **) calloc(sizeof(char *), k*num_queries);
  if (!all_ids) {
    goto malloc_failure;
  }
  all_scores = (double *) calloc(sizeof(double), k*num_queries);
  if (!all_scores) {
    goto malloc_failure;
  }
 skip_malloc:

  knearest_search->query_start = query_arena + (query_start*query_storage_size);
  knearest_search->num_queries = num_queries;
  knearest_search->query_fp_size = (num_bits+7)/8;
  knearest_search->query_storage_size = query_storage_size;
  
  knearest_search->k = k;
  knearest_search->search_state = 0;
  knearest_search->threshold = threshold;

  knearest_search->heaps = heaps;

  for (i=0; i<num_queries; i++) {
    heaps[i].indicies = all_indicies+(i*k);
    heaps[i].ids = all_ids+(i*k);
    heaps[i].scores = all_scores+(i*k);
  }
  knearest_search->num_targets_processed = 0;
  knearest_search->_all_ids = all_ids;
  knearest_search->_all_scores = all_scores;

  return CHEMFP_OK;


 malloc_failure:
  if (all_scores) free(all_scores);
  if (all_ids) free(all_ids);
  if (all_indicies) free(all_indicies);
  if (heaps) free(heaps);
  return CHEMFP_NO_MEM;
}

static char *new_string(const char *start, const char *end) {
  int n = end-start;
  char *s = malloc(n+1);
  if (s) {
    memcpy(s, start, n);
    s[n] = '\0';
  }
  return s;
}


int chemfp_fps_knearest_search_feed(
	chemfp_fps_knearest_search *knearest_search,
	int target_block_len, const char *target_block) {
  int k;
  double score, threshold;
  int num_added = 0;
  char *s;
  const char *line, *next_line, *end, *id_start, *id_end;
  const unsigned char *query_fp;
  chemfp_fps_heap *heap;
  int query_hex_size, query_fp_size;
  int i, err, retval;
  
  if (target_block_len == 0 || target_block[target_block_len-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;
  end = target_block+target_block_len;

  threshold = knearest_search->threshold;
  k = knearest_search->k;
  query_fp_size = knearest_search->query_fp_size;
  query_hex_size = query_fp_size * 2;

  line = target_block;
  while (line < end) {
    err = chemfp_fps_find_id(query_hex_size, line, &id_start, &id_end);
    if (err < 0) {
      retval = err;
      goto finish;
    }
    next_line = chemfp_to_next_line(id_end);
    query_fp = knearest_search->query_start;
    
    heap = knearest_search->heaps;
    for (i=0; i<knearest_search->num_queries; i++, query_fp += query_fp_size, heap++) {
      switch(heap->heap_state) {
      case ADD_TO_HEAP:
	score = chemfp_byte_hex_tanimoto(query_fp_size, query_fp, line);
	if (score >= threshold) {
	  heap->scores[heap->size] = score;
	  s = new_string(id_start, id_end);
	  if (!s) {
	    retval = CHEMFP_NO_MEM;
	    goto finish;
	  }
	  heap->ids[heap->size] = s;
	  heap->size++;
	}
	if (heap->size == k) {
	  chemfp_heapq_heapify(k, (void *)heap, (chemfp_heapq_lt) fps_heap_lt,
			       (chemfp_heapq_swap) fps_heap_swap);
	  heap->heap_state = REPLACE_IN_HEAP;
	}
	break;

      case REPLACE_IN_HEAP:
	score = chemfp_byte_hex_tanimoto(query_fp_size, query_fp, line);
	if (score > heap->scores[0]) {
	  heap->scores[0] = score;
	  free(heap->ids[0]);
	  s = new_string(id_start, id_end);
	  if (!s) {
	    retval = CHEMFP_NO_MEM;
	    goto finish;
	  }
	  heap->ids[0] = s;
	  chemfp_heapq_siftup(k, (void *) heap, 0,
			      (chemfp_heapq_lt) fps_heap_lt,
			      (chemfp_heapq_swap) fps_heap_swap);
	  if (heap->scores[0] == 1.0) {
	    heap->heap_state = MAXED_OUT_HEAP;
	  }
	}
	break;
      case MAXED_OUT_HEAP:
	continue;
      default:
	return -1; /* Not possible */
      }
    }
    line = next_line;
    num_added++;
  }
  retval = CHEMFP_OK;
 finish:
  knearest_search->num_targets_processed += num_added;
  return retval;
}

void chemfp_fps_knearest_search_free(chemfp_fps_knearest_search *knearest_search) {
  free(knearest_search->_all_scores);
  free(knearest_search->_all_ids);
  free(knearest_search->heaps);
}

void chemfp_fps_knearest_search_finish(chemfp_fps_knearest_search *knearest_search) {
  int i;
  chemfp_fps_heap *heap;
  if (knearest_search->search_state == 1) {
    return;
  }
  knearest_search->search_state = 1;
  for (i=0; i<knearest_search->num_queries; i++) {
    heap = knearest_search->heaps+i;
    if (heap->size < knearest_search->k) {
      chemfp_heapq_heapify(heap->size, (void *)heap,
                         (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
    }
    chemfp_heapq_heapsort(heap->size, (void *)heap,
			  (chemfp_heapq_lt) fps_heap_lt, (chemfp_heapq_swap) fps_heap_swap);
  }
}
