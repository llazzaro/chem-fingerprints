/* Functions for using the "fps" hex-based fingerprint file format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "heapq.h"
#include "chemfp.h"

enum {ADD_TO_HEAP, REPLACE_IN_HEAP, MAXED_OUT_HEAP};

/* Internal function to find the id field in an FPS line */
/* (Which means the fingerprint field is from line to *id_start-1 ) */
/* The line MUST match /^[0-9A-Fa-f]+\t[^\t\r\n]+/ */
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
  fp_field_len = (int) strspn(line, "0123456789abcdefABCDEF");
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
  case ' ': return CHEMFP_UNSUPPORTED_WHITESPACE;
  default: return CHEMFP_BAD_FINGERPRINT;
  }
  s++;

  /* You must pass in a newline-terminated string to this function.
     Therefore, this function will finish while inside the string.
     Note that I'm also checking for illegal whitespace here. */
  id_len = (int) strcspn(s, "\t\n\r");
  switch (s[id_len]) {
  case '\0': return CHEMFP_BAD_ID;
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

/* Return the number of fingerprints in the fps block which are greater
   than or equal to the specified threshold. */
int chemfp_fps_count_tanimoto_hits(
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
         query_index++, query_fp += query_storage_size) {
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
    /* There's no guarantee that the missing newline is on "stopped_at" */
    /* In the Python API there's no way to trigger this through normal code. */
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
         query_index++, query_fp += query_storage_size) {
      score = chemfp_byte_hex_tanimoto(fp_size, query_fp, line);
      if (score >= threshold) {
        current_cell->score = score;
        current_cell->query_index = query_index;
        current_cell->id_start = (int)(id_start - target_block);
        current_cell->id_end = (int)(id_end - target_block);
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
  *num_cells_processed = (int)(current_cell - cells);
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
  /* break ties on a first-come basis */
  return (heap->indices[i] > heap->indices[j]);
}

/* Swap two entries in the heap */
static void fps_heap_swap(chemfp_fps_heap *heap, int i, int j) {
  int idx = heap->indices[i];
  double score = heap->scores[i];
  char *id = heap->ids[i];

  heap->indices[i] = heap->indices[j];
  heap->scores[i] = heap->scores[j];
  heap->ids[i] = heap->ids[j];

  heap->indices[j] = idx;
  heap->scores[j] = score;
  heap->ids[j] = id;
}

/***************** new code */

int chemfp_fps_knearest_search_init(
        chemfp_fps_knearest_search *knearest_search,
        int num_bits, int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,
        int k, double threshold) {

  chemfp_fps_heap *heaps = NULL;
  int *all_indices = NULL;
  char **all_ids = NULL;
  double *all_scores = NULL;
  int i, num_queries;

  if (query_start >= query_end) {
    num_queries = 0;
    goto skip_malloc;
  } else {
    num_queries = query_end - query_start;
  }
  heaps = (chemfp_fps_heap *) calloc(num_queries, sizeof(chemfp_fps_heap));
  if (!heaps) {
    goto malloc_failure;
  }
  all_indices = (int *) calloc(k*num_queries, sizeof(int));
  if (!all_indices) {
    goto malloc_failure;
  }
  all_ids = (char **) calloc(k*num_queries, sizeof(char *));
  if (!all_ids) {
    goto malloc_failure;
  }
  all_scores = (double *) calloc(k*num_queries, sizeof(double));
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
    heaps[i].indices = all_indices+(i*k);
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
  if (all_indices) free(all_indices);
  if (heaps) free(heaps);
  return CHEMFP_NO_MEM;
}

static char *new_string(const char *start, const char *end) {
  size_t n = end-start;
  char *s = malloc(n+1);
  if (s) {
    memcpy(s, start, n);
    s[n] = '\0';
  }
  return s;
}


int chemfp_fps_knearest_tanimoto_search_feed(
        chemfp_fps_knearest_search *knearest_search,
        int target_block_len, const char *target_block) {
  int k;
  double score, threshold;
  int num_added = 0;
  char *s;
  const char *line, *next_line, *end, *id_start, *id_end;
  const unsigned char *query_fp;
  chemfp_fps_heap *heap;
  int query_hex_size, query_fp_size, query_storage_size;
  int i, err, retval;
  
  if (target_block_len == 0 || target_block[target_block_len-1] != '\n')
    return CHEMFP_MISSING_NEWLINE;
  end = target_block+target_block_len;

  threshold = knearest_search->threshold;
  k = knearest_search->k;
  query_fp_size = knearest_search->query_fp_size;
  query_hex_size = query_fp_size * 2;
  query_storage_size = knearest_search->query_storage_size;

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
    for (i=0; i<knearest_search->num_queries; i++, query_fp += query_storage_size, heap++) {
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
