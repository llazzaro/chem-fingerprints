#include <stdlib.h>
#include <stdio.h>

#include "chemfp.h"
#include "chemfp_internal.h"


chemfp_threshold_result *chemfp_alloc_threshold_results(int size) {
  /* Initializes all of the counts to 0 (and makes everything else nice too) */
  return (chemfp_threshold_result *) calloc(size, sizeof(chemfp_threshold_result));
}

void chemfp_free_results(int num_results, chemfp_threshold_result *results) {
  int i, j, num_hits;
  chemfp_hit_block *p, *q;
  for (i=0; i<num_results; i++) {
    num_hits = results[i].num_hits;
    if (num_hits > 1) {
      p = results[i].first;
      for (j=1+CHEMFP_HIT_BLOCK_SIZE; j<num_hits; j+=CHEMFP_HIT_BLOCK_SIZE) {
        q = p->next;
        free(p);
        p = q;
      }
      free(p);
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
  int i;
  chemfp_hit_block *block;
  if (num_hits == 0) {
    result->index0 = target_index;
    result->score0 = score;
    result->num_hits = 1;
    return 1;
  }
  
  i = (num_hits-1) % CHEMFP_HIT_BLOCK_SIZE;
  if (i == 0) {
    /* Need to allocate another block */
    block = (chemfp_hit_block *) malloc(sizeof(chemfp_hit_block));
    if (!block) {
      return 0;
    }
    if (num_hits == 1) {
      /* Add the first block; need to initialize start */
      result->first = block;
    } else {
      /* Otherwise, append to the end of the current end */
      result->last->next = block;
    }
    /* Link "end" to the new end block */
    result->last = block;

  } else {
    /* Don't need to allocate new space */
    block = result->last;
  }
  block->target_indices[i] = target_index;
  block->scores[i] = score;
  result->num_hits = num_hits+1;
  return 1;
}

int chemfp_threshold_result_get_hits(chemfp_threshold_result *result,
                                     chemfp_assign_hits_p add_callback, void *payload) {
  int num_hits = result->num_hits;
  int i, j, errval, index;
  chemfp_hit_block *block;

  if (num_hits) {
    add_callback(payload, 0, result->index0, result->score0);
    block = result->first;
    for (i=1; i<num_hits; i+=CHEMFP_HIT_BLOCK_SIZE) {
      for (j=0; j<CHEMFP_HIT_BLOCK_SIZE; j++) {
        index = i+j;
        if (index >= num_hits) {
          break;
        }
        errval = add_callback(payload, index, block->target_indices[j], block->scores[j]);
        if (errval) {
          return errval;
        }
      }
      block = block->next;
    }
  }
  return 0;
}
