#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"


chemfp_threshold_result *chemfp_alloc_threshold_results(int size) {
  /* Initializes all of the counts to 0 (and makes everything else nice too) */
  return (chemfp_threshold_result *) calloc(size, sizeof(chemfp_threshold_result));
}

void chemfp_free_results(int num_results, chemfp_threshold_result *results) {
  int i;
  for (i=0; i<num_results; i++) {
    if (results[i].num_hits) {
      free(results[i].scores);
    }
  }
  free(results);
}

int chemfp_get_num_hits(chemfp_threshold_result *result) {
  return result->num_hits;
}

int _chemfp_add_hit(chemfp_threshold_result *result,
                    int target_index, double score) {
  int num_hits = result->num_hits;
  int num_allocated = result->num_allocated;
  int *indices, *old_indices;
  double *scores;
  if (num_hits == num_allocated) {
    if (num_hits == 0) {
      num_allocated = 6;
      scores = (double *) malloc(num_allocated * (sizeof(int)+sizeof(double)));
      if (!scores) {
	return 0;
      }
      indices = (int *) (scores + num_allocated);
      result->num_allocated = num_allocated;
      result->indices = indices;
      result->scores = scores;
    } else {
      /* Grow by about 12% each time; this is the Python listobject resize strategy */
      num_allocated += (num_allocated >> 3) + (num_allocated < 9 ? 3 : 6);
      scores = (double *) realloc(result->scores, num_allocated * (sizeof(int)+sizeof(double)));
      if (!scores) {
	return 0;
      }
      /* Shift the indices to its new location */
      old_indices = (int *) (scores + num_hits);
      indices = (int *) (scores + num_allocated);
      memmove(indices, old_indices, num_hits*sizeof(int));
      result->num_allocated = num_allocated;
      result->indices = indices;
      result->scores = scores;
    }
  } else {
    indices = result->indices;
    scores = result->scores;
  }
  indices[num_hits] = target_index;
  scores[num_hits] = score;
  result->num_hits = num_hits+1;
  return 1;
}

int chemfp_threshold_result_get_hits(chemfp_threshold_result *result,
                                     chemfp_assign_hits_p add_callback, void *payload) {
  int num_hits = result->num_hits;
  int i, errval;

  if (num_hits) {
    for (i=0; i<num_hits; i++) {
      errval = add_callback(payload, i, result->indices[i], result->scores[i]);
      if (errval) {
	return errval;
      }
    }
  }
  return 0;
}
