#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "heapq.h"
#include "chemfp.h"

typedef struct {
  int *indicies;
  double *scores;
} IndexScoreData;

static int double_score_lt(IndexScoreData *data, int i, int j) {
  if (data->scores[i] < data->scores[j])
    return 1;
  if (data->scores[i] > data->scores[j])
    return 0;
  /* Sort in descending order by index. (XXX important or overkill?) */
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

enum scoring_directions {
  UP_OR_DOWN = 0,
  UP_ONLY, 
  DOWN_ONLY,
  FINISHED
};

typedef struct {
  int direction;
  int query_popcount;
  int max_popcount;
  int popcount;
  int up_popcount;
  int down_popcount;
  double score;
} PopcountSearchOrder;

static void init_search_order(PopcountSearchOrder *popcount_order, int query_popcount,
			      int max_popcount) {
  popcount_order->query_popcount = query_popcount;
  popcount_order->popcount = query_popcount;
  popcount_order->max_popcount = max_popcount;
  if (query_popcount <= 1) {
    popcount_order->direction = UP_ONLY;
    popcount_order->down_popcount = 0;
  } else {
    popcount_order->direction = UP_OR_DOWN;
    popcount_order->down_popcount = query_popcount-1;
  }
  popcount_order->up_popcount = query_popcount;
}

static void ordering_no_higher(PopcountSearchOrder *popcount_order) {
  switch (popcount_order->direction) {
  case UP_OR_DOWN:
    popcount_order->direction = DOWN_ONLY;
    break;
  case UP_ONLY:
    popcount_order->direction = FINISHED;
    break;
  default:
    break;
  }
}
static void ordering_no_lower(PopcountSearchOrder *popcount_order) {
  switch (popcount_order->direction) {
  case UP_OR_DOWN:
    popcount_order->direction = UP_ONLY;
    break;
  case DOWN_ONLY:
    popcount_order->direction = FINISHED;
    break;
  default:
    break;
  }
}


#define UP_SCORE(po) (((double)(po->query_popcount))/po->up_popcount)
#define DOWN_SCORE(po) (((double)(po->down_popcount))/po->query_popcount)

static int next_popcount(PopcountSearchOrder *popcount_order, double threshold) {
  double up_score, down_score;

  switch (popcount_order->direction) {
  case UP_OR_DOWN:
    up_score = UP_SCORE(popcount_order);
    down_score = DOWN_SCORE(popcount_order);
    if (up_score >= down_score) {
      popcount_order->popcount = (popcount_order->up_popcount)++;
      popcount_order->score = up_score;
      if (popcount_order->up_popcount > popcount_order->max_popcount) {
	popcount_order->direction = DOWN_ONLY;
      }
    } else {
      popcount_order->popcount = (popcount_order->down_popcount)--;
      popcount_order->score = down_score;
      if (popcount_order->down_popcount < 0) {
	popcount_order->direction = UP_ONLY;
      }
    }
    break;
   
  case UP_ONLY:
    popcount_order->score = UP_SCORE(popcount_order);
    popcount_order->popcount = (popcount_order->up_popcount)++;
    if (popcount_order->up_popcount > popcount_order->max_popcount) {
      popcount_order->direction = FINISHED;
    }
    break;
    
  case DOWN_ONLY:
    popcount_order->score = DOWN_SCORE(popcount_order);
    popcount_order->popcount = (popcount_order->down_popcount)--;
    if (popcount_order->down_popcount < 0) {
      popcount_order->direction = FINISHED;
    }
    break;

  default:
    return 0;
  }

  /* If the best possible score is under the threshold then we're done. */
  if (popcount_order->score < threshold) {
    popcount_order->direction = FINISHED;
    return 0;
  }
  return 1;

}

static int 
check_bounds(PopcountSearchOrder *popcount_order,
	     int *start, int *end, int target_start, int target_end) {
  if (*start > target_end) {
    ordering_no_higher(popcount_order);
    return 0;
  }
  if (*end < target_start) {
    ordering_no_lower(popcount_order);
    return 0;
  };

  if (*start < target_start) {
    ordering_no_higher(popcount_order);
    *start = target_start;
  }
  if (*end > target_end) {
    ordering_no_lower(popcount_order);
    *end = target_end;
  }
  return 1;
}

/* count code */
int chemfp_count_tanimoto_arena(
	/* Count all matches within the given threshold */
	double threshold,

	/* Number of bits in the fingerprint */
	int num_bits,

	/* Query arena, start and end indicies */
	int query_storage_size,
	const unsigned char *query_arena, int query_start, int query_end,

	/* Target arena, start and end indicies */
	int target_storage_size,
	const unsigned char *target_arena, int target_start, int target_end,

	/* Target popcount distribution information */
	int *target_popcount_indicies,

	/* Results go into these arrays  */
	int *result_counts
				   ) {

  int query_index, target_index;
  const unsigned char *query_fp, *target_fp;
  int start, end;
  int count;
  int fp_size = num_bits / 8;
  double score, popcount_sum;
  int query_popcount, start_target_popcount, end_target_popcount;
  int target_popcount;
  int intersect_popcount;
  
  if (query_start >= query_end) {
    /* No queries */
    return 0;
  }
  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  if (threshold < 1.0/num_bits) {
    threshold = 0.0;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    for (query_index = query_start; query_index < query_end; query_index++) {
      /* No possible targets */
      *result_counts++ = 0;
    }
    return query_index-query_start;
  }
  if (target_popcount_indicies == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
    query_fp = query_arena + (query_start * query_storage_size);
    for (query_index = query_start; query_index < query_end;
	 query_index++, query_fp += query_storage_size) {
      target_fp = target_arena + (target_start * target_storage_size);
      /* Handle the popcount(query) == 0 special case? */
      count = 0;

      for (target_index = target_start; target_index < target_end;
	   target_index++, target_fp += query_storage_size) {
	score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
	if (score >= threshold) {
	  count++;
	}
      }
      *result_counts++ = count;
    }
    return query_index-query_start;
  }
  
  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the ordering because it's supposed to find everything */

  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {
    
    query_popcount = chemfp_byte_popcount(fp_size, query_fp);
    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
	*result_counts++ = (target_end - target_start);
      } else {
	result_counts++;
      }
      continue;
    }
    /* Figure out which fingerprints to search */
    if (threshold == 0.0) {
      start_target_popcount = 0;
      end_target_popcount = num_bits;
    } else {
      start_target_popcount = query_popcount * threshold;
      end_target_popcount = ceil(query_popcount / threshold);
      if (end_target_popcount > num_bits) {
	end_target_popcount = num_bits;
      }
    }
    count = 0;
    for (target_popcount=start_target_popcount; target_popcount<=end_target_popcount;
	 target_popcount++) {
      start = target_popcount_indicies[target_popcount];
      end = target_popcount_indicies[target_popcount+1];
      if (start < target_start) {
	start = target_start;
      }
      if (end > target_end) {
	end = target_end;
      }

      target_fp = target_arena + (start * target_storage_size);
      popcount_sum = query_popcount + target_popcount;
      for (target_index = start; target_index < end;
	   target_index++, target_fp += target_storage_size) {
	intersect_popcount = chemfp_byte_intersect_popcount(fp_size, query_fp, target_fp);
	score = intersect_popcount / (popcount_sum - intersect_popcount);
	if (score >= threshold) {
	  count++;
	}
      }
    }
    *result_counts++ = count;
  } /* went through each of the queries */
  return query_index-query_start;
}


int chemfp_threshold_tanimoto_arena(
	/* Within the given threshold */
	double threshold,

	/* Number of bits in the fingerprint */
	int num_bits,

	/* Query arena, start and end indicies */
	int query_storage_size, const unsigned char *query_arena,
	int query_start, int query_end,

	/* Target arena, start and end indicies */
	int target_storage_size, const unsigned char *target_arena,
	int target_start, int target_end,

	/* Target popcount distribution information */
	/*  (must have at least num_bits+1 elements) */
	int *target_popcount_indicies,

	/* Results go into these arrays  */
	int *result_offsets,
	int num_cells,
	int *result_indicies,
	double *result_scores
				    ) {

  int query_index, target_index;
  const unsigned char *query_fp, *target_fp;
  int start, end;
  int count;
  int fp_size = (num_bits+7) / 8;
  double score;
  int query_popcount, start_target_popcount, end_target_popcount;
  int target_popcount;
  int intersect_popcount, popcount_sum;
  int result_offset = *result_offsets++;
  
  if (query_start >= query_end) {
    /* No queries */
    return 0;
  }

  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  if (threshold < 1.0/num_bits) {
    threshold = 0.0;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    for (query_index = query_start; query_index < query_end; query_index++) {
      /* No possible targets */
      *result_offsets++ = result_offset;
    }
    return query_index-query_start;
  }
  if (target_popcount_indicies == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
    query_fp = query_arena + (query_start * query_storage_size);
    for (query_index = query_start; query_index < query_end;
	 query_index++, query_fp += query_storage_size) {
      target_fp = target_arena + (target_start * target_storage_size);
      /* Handle the popcount(query) == 0 special case? */
      count = 0;
      for (target_index = target_start; target_index < target_end;
	   target_index++, target_fp += query_storage_size) {
	score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
	if (score >= threshold) {
	  *result_indicies++ = target_index;
	  *result_scores++ = score;
	  count++;
	}
      }
      result_offset += count;
      num_cells -= count;
      *result_offsets++ = result_offset;
    }
    return query_index-query_start;
  }
  
  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the ordering because it's supposed to find everything */

  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {

    if (num_cells < (target_end - target_start)) {
      break;
    }
    
    query_popcount = chemfp_byte_popcount(fp_size, query_fp);
    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
	for (target_index = target_start; target_index < target_end; target_index++) {
	  *result_indicies++ = target_index;
	  *result_scores++ = 0.0;
	}
	count = (target_index - target_start);
	result_offset += count;
	num_cells -= count;
	*result_offsets++ = result_offset;
      }
      continue;
    }
    /* Figure out which fingerprints to search */
    if (threshold == 0.0) {
      start_target_popcount = 0;
      end_target_popcount = num_bits;
    } else {
      start_target_popcount = query_popcount * threshold;
      end_target_popcount = ceil(query_popcount / threshold);
      if (end_target_popcount > num_bits) {
	end_target_popcount = num_bits;
      }
    }

    count = 0;
    for (target_popcount=start_target_popcount; target_popcount<=end_target_popcount;
	 target_popcount++) {
      start = target_popcount_indicies[target_popcount];
      end = target_popcount_indicies[target_popcount+1];
      if (start < target_start) {
	start = target_start;
      }
      if (end > target_end) {
	end = target_end;
      }

      target_fp = target_arena + (start * target_storage_size);
      popcount_sum = query_popcount + target_popcount;
      for (target_index = start; target_index < end;
	   target_index++, target_fp += target_storage_size) {
	intersect_popcount = chemfp_byte_intersect_popcount(fp_size, query_fp, target_fp);
	score = ((double) intersect_popcount) / (popcount_sum - intersect_popcount);
	if (score >= threshold) {
	  *result_indicies++ = target_index;
	  *result_scores++ = score;
	  count++;
	}
      }
    }
    result_offset += count;
    num_cells -= count;
    *result_offsets++ = result_offset;
  } /* went through each of the queries */
  return query_index-query_start;
}





/**** k-nearest code ****/

static int 
knearest_tanimoto_arena_no_popcounts(
	/* Find the 'k' nearest items */
	int k,
	/* Within the given threshold */
	double threshold,

	/* Fingerprint size in bits */
	int num_bits,

	/* Query arena, start and end indicies */
	int query_storage_size, const unsigned char *query_arena,
	int query_start, int query_end,

	/* Target arena, start and end indicies */
	int target_storage_size, const unsigned char *target_arena,
	int target_start, int target_end,

	/* Results go into these arrays  */
	int *result_offsets,
	int num_cells,
	int *result_indicies,
	double *result_scores
				   ) {
  int query_index, target_index;
  int fp_size = (num_bits+7)/8;
  const unsigned char *query_fp, *target_fp;
  double query_threshold, score;
  int heap_size;
  IndexScoreData heap;
  int result_offset=*result_offsets++;

  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {

    if (num_cells < k) {
      /* Not enough space to store everything, so stop here. */
      return query_index-query_start;
    }

    query_threshold = threshold;
    heap_size = 0;
    heap.indicies = result_indicies;
    heap.scores = result_scores;
    
    target_fp = target_arena + (target_start * query_storage_size);
    target_index = target_start;

    for (; target_index < target_end;
	 target_index++, target_fp += target_storage_size) {
      score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
      if (score >= query_threshold) {
	heap.indicies[heap_size] = target_index;
	heap.scores[heap_size] = score;
	heap_size++;
	if (heap_size == k) {
	  chemfp_heapq_heapify(heap_size, &heap,  (chemfp_heapq_lt) double_score_lt,
			       (chemfp_heapq_swap) double_score_swap);
	  query_threshold = heap.scores[0];
	  /* Since we leave the loop early, I need to advance the pointers */
	  target_index++;
	  target_fp += target_storage_size;
	  break;
	}
      }
    }
    /* Either we've reached the end of the fingerprints or the heap is full */
    if (heap_size == k) {
      /* Continue scanning through the fingerprints */
      for (; target_index < target_end;
	   target_index++, target_fp += target_storage_size) {
	score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);

	/* We need to be strictly *better* than what's in the heap */
	if (score > query_threshold) {
	  heap.indicies[0] = target_index;
	  heap.scores[0] = score;
	  chemfp_heapq_siftup(heap_size, &heap, 0, (chemfp_heapq_lt) double_score_lt,
			      (chemfp_heapq_swap) double_score_swap);
	  query_threshold = heap.scores[0];
	} /* heapreplaced the old smallest item with the new item */
      }
      /* End of the fingerprint scan */
    } else {
      /* The heap isn't full, so we haven't yet heapified it. */
      chemfp_heapq_heapify(heap_size, &heap,  (chemfp_heapq_lt) double_score_lt,
			   (chemfp_heapq_swap) double_score_swap);
    }

    /* Sort the elements */
    chemfp_heapq_heapsort(heap_size, &heap, (chemfp_heapq_lt) double_score_lt,
			  (chemfp_heapq_swap) double_score_swap);
    
    /* Pass back the query results */
    result_offset += heap_size;

    *result_offsets++ = result_offset;
    memcpy(result_indicies, heap.indicies, heap_size * sizeof(int));
    memcpy(result_scores, heap.scores, heap_size * sizeof(double));

    /* Move the pointers so I can report the next results */
    result_indicies += heap_size;
    result_scores += heap_size;
    num_cells -= heap_size;
  } /* Loop through the queries */

  return query_index-query_start;
}


int chemfp_knearest_tanimoto_arena(
	/* Find the 'k' nearest items */
	int k,
	/* Within the given threshold */
	double threshold,

	/* Size of the fingerprints and size of the storage block */
	int num_bits,

	/* Query arena, start and end indicies */
	int query_storage_size, const unsigned char *query_arena,
	int query_start, int query_end,

	/* Target arena, start and end indicies */
	int target_storage_size, const unsigned char *target_arena,
	int target_start, int target_end,

	/* Target popcount distribution information */
	int *target_popcount_indicies,

	/* Results go into these arrays  */
	int *result_offsets,
	int num_cells,
	int *result_indicies,
	double *result_scores
				   ) {

  int heap_size, fp_size;
  int query_popcount, target_popcount, intersect_popcount;
  double score, best_possible_score, popcount_sum, query_threshold;
  const unsigned char *query_fp, *target_fp;
  int query_index, target_index;
  int start, end;
  PopcountSearchOrder popcount_order;
  IndexScoreData heap;
  int result_offset;

  /* This is C. We don't check for illegal input values. */

  if (query_start >= query_end) {
    return 0;
  }

  /* k == 0 is a valid input, and of course the result is no matches */
  if (k == 0) {
    result_offset = *result_offsets++;
    for (query_index = query_start; query_index < query_end; query_index++) {
      *result_offsets++ = result_offset;
    }
    return query_index-query_start;
  }
  fp_size = (num_bits+7)/8;

  if (target_popcount_indicies == NULL) {
    /* precomputed targets aren't available. Use the slower algorithm. */
    return knearest_tanimoto_arena_no_popcounts(
	k, threshold, fp_size,
	query_storage_size, query_arena, query_start, query_end,
	target_storage_size, target_arena, target_start, target_end,
	result_offsets, num_cells, result_indicies, result_scores);
  }
  result_offset = *result_offsets++;

  /* Loop through the query fingerprints */
  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {
    if (num_cells < k) {
      /* Not enough space to store everything, so stop here. */
      return query_index-query_start;
    }

    query_threshold = threshold;
    query_popcount = chemfp_byte_popcount(fp_size, query_fp);

    if (query_popcount == 0) {
      /* By definition this will never return hits. Even if threshold == 0.0. */
      /* (I considered returning the first k hits, but that's chemically meaninless.) */
      /* XXX change this. Make it returns the first k hits */
      *result_offsets++ = result_offset;
      continue;
    }

    /* Search the bins using the ordering from Swamidass and Baldi.*/
    init_search_order(&popcount_order, query_popcount, num_bits);

    heap_size = 0;
    heap.indicies = result_indicies;
    heap.scores = result_scores;

    /* Look through the sections of the arena in optimal popcount order */
    while (next_popcount(&popcount_order, query_threshold)) {
      target_popcount = popcount_order.popcount;
      best_possible_score = popcount_order.score;

      /* If we can't beat the query threshold then we're done with the targets */
      if (best_possible_score < query_threshold) {
	break;
      }

      /* Scan through the targets which have the given popcount */
      start = target_popcount_indicies[target_popcount];
      end = target_popcount_indicies[target_popcount+1];
      
      if (!check_bounds(&popcount_order, &start, &end, target_start, target_end)) {
	continue;
      }

      /* Iterate over the target fingerprints */
      target_fp = target_arena + start*target_storage_size;
      popcount_sum = (double)(query_popcount + target_popcount);

      target_index = start;

      /* There are fewer than 'k' elements in the heap*/
      if (heap_size < k) {
	for (; target_index<end; target_index++, target_fp += target_storage_size) {
	  intersect_popcount = chemfp_byte_intersect_popcount(fp_size, query_fp, target_fp);
	  score = intersect_popcount / (popcount_sum - intersect_popcount);

	  /* The heap isn't full; only check if we're at or above the query threshold */
	  if (score >= query_threshold) {
	    heap.indicies[heap_size] = target_index;
	    heap.scores[heap_size] = score;
	    heap_size++;
	    if (heap_size == k) {
	      chemfp_heapq_heapify(heap_size, &heap,  (chemfp_heapq_lt) double_score_lt,
				   (chemfp_heapq_swap) double_score_swap);
	      query_threshold = heap.scores[0];
	      /* We're going to jump to the "heap is full" section */
	      /* Since we leave the loop early, I need to advance the pointers */
	      target_index++;
	      target_fp += target_storage_size;
	      goto heap_replace;
	    }
	  } /* Added to heap */
	} /* Went through target fingerprints */

	/* If we're here then the heap did not fill up. Try the next popcount */
	continue;
      }

    heap_replace:
      /* We only get here if the heap contains k element */

      /* Earlier we tested for "best_possible_score<query_threshold". */
      /* The test to replace an element in the heap is more stringent. */
      if (query_threshold >= best_possible_score) {
	/* Can't do better. Might as well give up. */
	break;
      }

      /* Scan through the target fingerprints; can we improve over the threshold? */
      for (; target_index<end; target_index++, target_fp += target_storage_size) {
	intersect_popcount = chemfp_byte_intersect_popcount(fp_size, query_fp, target_fp);
	score = intersect_popcount / (popcount_sum - intersect_popcount);

	/* We need to be strictly *better* than what's in the heap */
	if (score > query_threshold) {
	  heap.indicies[0] = target_index;
	  heap.scores[0] = score;
	  chemfp_heapq_siftup(heap_size, &heap, 0, (chemfp_heapq_lt) double_score_lt,
			      (chemfp_heapq_swap) double_score_swap);
	  query_threshold = heap.scores[0];
	  if (query_threshold >= best_possible_score) {
	    /* we can't do any better in this section (or in later ones) */
	    break;
	  }
	} /* heapreplaced the old smallest item with the new item */
      } /* looped over fingerprints */
    } /* Went through all the popcount regions */

    /* We have scanned all the fingerprints. Is the heap full? */

    if (heap_size < k) {
      /* Not full, so need to heapify it. */
      chemfp_heapq_heapify(heap_size, &heap,  (chemfp_heapq_lt) double_score_lt,
			   (chemfp_heapq_swap) double_score_swap);
    }
    /* Sort the elements */
    chemfp_heapq_heapsort(heap_size, &heap, (chemfp_heapq_lt) double_score_lt,
			  (chemfp_heapq_swap) double_score_swap);
    
    /* Pass back the query results */
    result_offset += heap_size;
    *result_offsets++ = result_offset;
    if (heap_size) {
      result_indicies += heap_size;
      result_scores += heap_size;
      num_cells -= heap_size;
    }

  } /* looped over all queries */

  return query_index-query_start;
}

static int compare_by_popcount(const void *left_p, const void *right_p) {
  const ChemFPOrderedPopcount *left = (ChemFPOrderedPopcount *) left_p;
  const ChemFPOrderedPopcount *right = (ChemFPOrderedPopcount *) right_p;
  if (left->popcount < right->popcount) {
    return -1;
  }
  if (left->popcount > right->popcount) {
    return 1;
  }
  if (left->index < right->index) {
    return -1;
  }
  if (left->index > right->index) {
    return 1;
  }
  return 0;
}

int
chemfp_reorder_by_popcount(
	int num_bits,
	int storage_size, const unsigned char *arena, int start, int end,
	unsigned char *new_arena, ChemFPOrderedPopcount *ordering,
	int *popcount_indicies) {

  int num_fingerprints, popcount;
  int fp_size = (num_bits+7)/8;
  int fp_index, i;
  const unsigned char *fp;

  if (start >= end) {
    return 0;
  }

  num_fingerprints = end - start;

  fp = arena + (start * storage_size);
  for (fp_index = start; fp_index < end; fp_index++, fp += storage_size) {
    popcount = chemfp_byte_popcount(fp_size, fp);
    ordering[fp_index].popcount = popcount;
    ordering[fp_index].index = fp_index;
  }
  qsort(ordering, num_fingerprints, sizeof(ChemFPOrderedPopcount), compare_by_popcount);


  /* Build the new arena based on the values in the old arena */
  for (i=0; i<num_fingerprints; i++) {
    memcpy(new_arena, arena+(ordering[i].index * storage_size), storage_size);
    new_arena += storage_size;
  }

  /* Create the popcount indicies */
  if (popcount_indicies != NULL) {
    /* Since we've sorted by popcount, this is easy */
    popcount = 0;
    *popcount_indicies++ = 0;
    for (i=0; i<num_fingerprints; i++) {
      while (popcount < ordering[i].popcount) {
	*popcount_indicies++ = i;
	popcount++;
	if (popcount == num_bits) {
	  /* We are at or above the limit. We can stop now. */
	  i = num_fingerprints+1;
	  break;
	  /* Note: with corrupted data it is possible
	     that ->popcount can be > num_bits. This is
	     undefined behavior. I get to do what I want.
	     I decided to treat them as having "max_popcount" bits.
	     After all, I don't want corrupt data to crash the
	     system, and no one is going to validate the input
	     fingerprints for correctness each time.  */
	}
      }
    }
    /* Finish up the high end */
    while (popcount <= num_bits) {
      *popcount_indicies++ = i;
      popcount++;
    }
  }
  return num_fingerprints;
}
