#include <limits.h>
#include "heapq.h"
#include "chemfp.h"

typedef struct {
  int *indicies;
  double *scores;
} DoubleScoreData;

// XXX Check that heapq preserves order
int double_score_lt(DoubleScoreData *data, int i, int j) {
  if (data->scores[i] < data->scores[j])
	return 1;
  if (data->scores[i] > data->scores[j])
	return 0;
  // Sort in descending order by index. (XXX immportant or overkill?)
  return (data->indicies[i] >= data->indicies[j]);
}
void double_score_swap(DoubleScoreData *data, int i, int j) {
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
  DoubleScoreData heap;
  heap.indicies = indicies;
  heap.scores = scores;

  for (fp_index=0; num_added<n && fp_index<num_targets; fp_index++) {
	score = chemfp_byte_tanimoto(len, query_fp, target_block);
	if (score >= threshold) {
	  indicies[num_added] = fp_index;
	  scores[num_added] = score;
	  target_block += storage_len;
	  num_added++;
	}
  }
  heapq_heapify(num_added, &heap,
				(heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
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
		heapq_siftup(n, &heap, 0,
					 (heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
		threshold = scores[0]; // Omitting this is hard to test
	  }
	  target_block += storage_len;
	  fp_index++;
	}
  }
  heapq_heapsort(n, &heap,
				 (heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
  return n;
}
