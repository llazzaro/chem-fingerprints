#include <limits.h>
#include "heapq.h"
#include "chemfp.h"

typedef struct {
  int *indicies;
  double *scores;
} DoubleScoreData;

// XXX Check that heapq preserves order
int double_score_lt(DoubleScoreData *data, int i, int j) {
  return data->scores[i] < data->scores[j];
}
void double_score_swap(DoubleScoreData *data, int i, int j) {
  int tmp_index = data->indicies[i];
  double tmp_score = data->scores[i];
  data->indicies[i] = data->indicies[j];
  data->scores[i] = data->scores[j];
  data->indicies[j] = tmp_index;
  data->scores[j] = tmp_score;
}

int chemfp_nlargest_tanimoto(
		int len, unsigned char *query_fp,
		int num_fps, unsigned char *fps, int offset, int storage_len,
		unsigned int n, int *indicies, double *scores) {
  int i, fp_index;

  if (len < 1) {
	return -1;
  }
  if (storage_len == -1) {
	storage_len = len;
  } else if (len > storage_len) {
	return -1;
  }
  if (num_fps < 0) {
	return -1;
  }
  if (((long long) num_fps * storage_len) > INT_MAX) {
	return -1;
  }
  fps += offset;

  for (i=0; i<n; i++) {
	indicies[i] = -1;
	scores[i] = -1.0;
  }
  /* Preamble done. Let's get to work. */
  DoubleScoreData heap;
  heap.indicies = indicies;
  heap.scores = scores;

  for (fp_index=0; fp_index<n && fp_index<num_fps; fp_index++) {
	indicies[fp_index] = fp_index;
	scores[fp_index] = chemfp_byte_tanimoto(len, query_fp, fps);
	fps += storage_len;
  }
  heapq_heapify(fp_index, &heap,
				(heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
  if (fp_index < n) {
	n = fp_index;
  } else {
	double lowest_score, score;
	lowest_score = scores[0];
	while (fp_index < num_fps) {
	  score = chemfp_byte_tanimoto(len, query_fp, fps);
	  if (lowest_score < score) {
		scores[0] = score;
		indicies[0] = fp_index;
		heapq_siftup(n, &heap, 0,
					 (heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
	  }
	  fps += storage_len;
	  fp_index++;
	}
  }
  heapq_heapsort(n, &heap,
				 (heapq_lt) double_score_lt, (heapq_swap) double_score_swap);
  return n;
}
