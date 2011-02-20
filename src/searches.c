#include <limits.h>
#include <ctype.h>
#include "heapq.h"
#include "chemfp.h"

typedef struct {
  int *indicies;
  double *scores;
} IndexScoreData;

// XXX Check that heapq preserves order
static int double_score_lt(IndexScoreData *data, int i, int j) {
  if (data->scores[i] < data->scores[j])
    return 1;
  if (data->scores[i] > data->scores[j])
    return 0;
  // Sort in descending order by index. (XXX immportant or overkill?)
  return (data->indicies[i] >= data->indicies[j]);
}
static void double_score_swap(IndexScoreData *data, int i, int j) {
  int tmp_index = data->indicies[i];
  double tmp_score = data->scores[i];
  data->indicies[i] = data->indicies[j];
  data->scores[i] = data->scores[j];
  data->indicies[j] = tmp_index;
  data->scores[j] = tmp_score;
}

int chemfp_nlargest_tanimoto_block(
        int n,
        int len, unsigned char *query_fp,
        int num_targets, unsigned char *target_block, int offset, int storage_len,
        double threshold,
        int *indicies, double *scores) {
  int fp_index;
  int num_added = 0;
  double score;

  if (n < 0 || len < 1) {
    return -1;
  }
  if (storage_len == -1) {
    storage_len = len;
  } else if (len > storage_len) {
    return -1;
  }
  if (num_targets < 0) {
    return -1;
  }
  if (((long long) num_targets * storage_len) > INT_MAX) {
    return -1;
  }
  target_block += offset;

  /* Preamble done. Let's get to work. */
  IndexScoreData heap;
  heap.indicies = indicies;
  heap.scores = scores;

  for (fp_index=0; num_added<n && fp_index<num_targets; fp_index++) {
    score = chemfp_byte_tanimoto(len, query_fp, target_block);
    if (score >= threshold) {
      indicies[num_added] = fp_index;
      scores[num_added] = score;
      num_added++;
    }
    target_block += storage_len;
  }
  chemfp_heapq_heapify(num_added, &heap, (chemfp_heapq_lt) double_score_lt,
                       (chemfp_heapq_swap) double_score_swap);
                
  if (num_added < n) {
    /* Stopped because there are no more targets */
    n = num_added;
  } else {
    /* Process the rest of the targets */
    /* Reset the threshold to the smallest value in the heap */
    threshold = scores[0];
    while (fp_index < num_targets) {
      score = chemfp_byte_tanimoto(len, query_fp, target_block);
      if (threshold < score) {
        scores[0] = score;
        indicies[0] = fp_index;
        chemfp_heapq_siftup(n, &heap, 0, (chemfp_heapq_lt) double_score_lt,
                            (chemfp_heapq_swap) double_score_swap);
        threshold = scores[0]; // Omitting this is hard to test
      }
      target_block += storage_len;
      fp_index++;
    }
  }
  chemfp_heapq_heapsort(n, &heap, (chemfp_heapq_lt) double_score_lt,
                        (chemfp_heapq_swap) double_score_swap);
                 
  return n;
}


static int hex_readline(int len, unsigned char *hex_query_fp,
                        unsigned char *target_block, int *offset_p,
                        double *score_p, unsigned char **start_id_p, int *id_len_p,
                        int *lineno) {

  int offset = *offset_p;

  // Must start with a hex fingerprint
  *score_p = chemfp_hex_tanimoto(len, hex_query_fp, target_block+offset);
  if (*score_p == -1.0) {
    // target was not a hex or wasn't long enough
    return -4;
  }
  // Go to the character after the hex string
  offset += len;
  // This must end inside of target_block because it terminates with an '\n'
  for (; isspace(target_block[offset]); offset++)
    ;
  // Either the end of the line or the start of the id block
  if (target_block[offset] == '\n') {
    return -4; // missing id
  }
  // Start of the id
  *start_id_p = target_block+offset;
  for (; !isspace(target_block[offset]); offset++)
    ;
  // End of the id
  *id_len_p = target_block+offset-*start_id_p;
  // Skip to end of line if not already there
  for (;target_block[offset] != '\n'; offset++)
    ;
  // Skip over the newline character; ready for the next line
  offset++;

  *offset_p = offset;
  (*lineno)++;
  return 0;
}

typedef struct {
  double *scores;
  unsigned char **start_ids;
  int *id_lens;
} HexScoreData;

static int hex_score_lt(HexScoreData *data, int i, int j) {
  if (data->scores[i] < data->scores[j])
    return 1;
  if (data->scores[i] > data->scores[j])
    return 0;
  return i < j; // XXX right?
}
static void hex_score_swap(HexScoreData *data, int i, int j) {
  double tmp_score = data->scores[i];
  unsigned char *tmp_start_id = data->start_ids[i];
  int tmp_id_len = data->id_lens[i];

  data->scores[i] = data->scores[j];
  data->start_ids[i] = data->start_ids[j];
  data->id_lens[i] = data->id_lens[j];

  data->scores[j] = tmp_score;
  data->start_ids[j] = tmp_start_id;
  data->id_lens[j] = tmp_id_len;
}

int chemfp_hex_tanimoto_block(
        int n,
        int len, unsigned char *hex_query_fp,
        int target_len, unsigned char *target_block,
        double threshold,
        double *scores, unsigned char **start_ids, int *id_lens, int *lineno_p) {
  double line_score;
  int line_id_len;
  unsigned char *line_start_id;


  HexScoreData heap;
  int err;
  heap.scores = scores;
  heap.start_ids = start_ids;
  heap.id_lens = id_lens;

  if (n < 0 || len < 1 || target_len < 0)
    return -1;
  if (target_block[target_len-1] != '\n')
    return -2;
  if (!chemfp_hex_isvalid(len, hex_query_fp))
    return -3;

  int offset = 0;
  int num_added = 0;
  while (num_added<n && offset<target_len) {
    // I know there's at least one to add
    err = hex_readline(len, hex_query_fp, target_block, &offset,
                       &line_score, &line_start_id, &line_id_len, lineno_p);
    if (err < 0)
      return err;
    // 
    if (line_score >= threshold) {
      heap.scores[num_added] = line_score;
      heap.start_ids[num_added] = line_start_id;
      id_lens[num_added] = line_id_len;
      num_added++;
    }
  }
  chemfp_heapq_heapify(num_added, &heap, (chemfp_heapq_lt) hex_score_lt,
                       (chemfp_heapq_swap) hex_score_swap);
                
  if (num_added < n) {
    /* Stopped because there are no more targets */
    n = num_added;
  } else {
    threshold = scores[0];
    while (offset < target_len) {
      err = hex_readline(len, hex_query_fp, target_block, &offset,
                         &line_score, &line_start_id, &line_id_len, lineno_p);
      if (err < 0)
        return err;
      if (threshold < line_score) {
        scores[0] = line_score;
        start_ids[0] = line_start_id;
        id_lens[0] = line_id_len;
        chemfp_heapq_siftup(n, &heap, 0, (chemfp_heapq_lt) hex_score_lt,
                            (chemfp_heapq_swap) hex_score_swap);
                     
        threshold = scores[0]; // Omitting this is hard to test
      }
    }
  }
  chemfp_heapq_heapsort(n, &heap, (chemfp_heapq_lt) hex_score_lt,
                        (chemfp_heapq_swap) hex_score_swap);
  return n;
}

/* Count the number of byte fingerprints which, when intersected with the
   query, have at least min_overlap bits in common */

int chemfp_byte_intersect_popcount_count(
        int len, unsigned char *query_fp,
        int num_targets, unsigned char *target_block, int offset, int storage_len,
		int min_overlap) {
  int count = 0;
  int fp_index;
  int popcount;
  if (len < 1) {
	return -1;
  }
  if (storage_len == -1) {
    storage_len = len;
  } else if (len > storage_len) {
    return -1;
  }
  if (num_targets < 0) {
    return -1;
  }
  if (((long long) num_targets * storage_len) > INT_MAX) {
    return -1;
  }
  target_block += offset;

  /* Let's count! */
  for (fp_index=0; fp_index < num_targets; fp_index++) {
	popcount = chemfp_byte_intersect_popcount(len, query_fp, target_block);
	if (popcount >= min_overlap) {
	  count += 1;
	}
	target_block += storage_len;
  }
  return count;
}
