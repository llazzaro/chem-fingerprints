#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdio.h> // XXX remove
#include <stdlib.h>
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
} PopcountOrdering;

static void init_popcount_ordering(PopcountOrdering *popcount_ordering, int query_popcount,
				   int max_popcount) {
  //  printf("init_popcount_ordering(%d, %d)\n", query_popcount, max_popcount);
  popcount_ordering->query_popcount = query_popcount;
  popcount_ordering->popcount = query_popcount;
  popcount_ordering->max_popcount = max_popcount;
  if (query_popcount <= 1) {
    popcount_ordering->direction = UP_ONLY;
    popcount_ordering->down_popcount = 0;
  } else {
    popcount_ordering->direction = UP_OR_DOWN;
    popcount_ordering->down_popcount = query_popcount-1;
  }
  popcount_ordering->up_popcount = query_popcount;
}
static void ordering_no_higher(PopcountOrdering *popcount_ordering) {
  switch (popcount_ordering->direction) {
  case UP_OR_DOWN:
    popcount_ordering->direction = DOWN_ONLY;
    break;
  case UP_ONLY:
    popcount_ordering->direction = FINISHED;
    break;
  default:
    break;
  }
}
static void ordering_no_lower(PopcountOrdering *popcount_ordering) {
  switch (popcount_ordering->direction) {
  case UP_OR_DOWN:
    popcount_ordering->direction = UP_ONLY;
    break;
  case DOWN_ONLY:
    popcount_ordering->direction = FINISHED;
    break;
  default:
    break;
  }
}


#define UP_SCORE(po) (((double)(po->query_popcount))/po->up_popcount)
#define DOWN_SCORE(po) (((double)(po->down_popcount))/po->query_popcount)

static int next_popcount(PopcountOrdering *popcount_ordering, double threshold) {
  double up_score, down_score;
  //  printf("Starting with %d %d up %d down %d\n", popcount_ordering->direction,
  //	 popcount_ordering->popcount,
  //	 popcount_ordering->up_popcount, popcount_ordering->down_popcount);
  switch (popcount_ordering->direction) {
  case UP_OR_DOWN:
    up_score = UP_SCORE(popcount_ordering);
    down_score = DOWN_SCORE(popcount_ordering);
    //printf("up %f down %f\n", up_score, down_score);
    if (up_score >= down_score) {
      popcount_ordering->popcount = (popcount_ordering->up_popcount)++;
      popcount_ordering->score = up_score;
      if (popcount_ordering->up_popcount > popcount_ordering->max_popcount) {
	popcount_ordering->direction = DOWN_ONLY;
      }
    } else {
      popcount_ordering->popcount = (popcount_ordering->down_popcount)--;
      popcount_ordering->score = down_score;
      if (popcount_ordering->down_popcount < 0) {
	popcount_ordering->direction = UP_ONLY;
      }
    }
    break;
   
  case UP_ONLY:
    popcount_ordering->score = UP_SCORE(popcount_ordering);
    popcount_ordering->popcount = (popcount_ordering->up_popcount)++;
    if (popcount_ordering->up_popcount > popcount_ordering->max_popcount) {
      popcount_ordering->direction = FINISHED;
    }
    break;
    
  case DOWN_ONLY:
    popcount_ordering->score = DOWN_SCORE(popcount_ordering);
    popcount_ordering->popcount = (popcount_ordering->down_popcount)--;
    if (popcount_ordering->down_popcount < 0) {
      popcount_ordering->direction = FINISHED;
    }
    break;

  default:
    return 0;
  }

  //  printf("count and score %d %f (%f)\n",
  // popcount_ordering->popcount, popcount_ordering->score, threshold);

  /* If the best possible score is under the threshold then we're done. */
  if (popcount_ordering->score < threshold) {
    popcount_ordering->direction = FINISHED;
    return 0;
  }
  return 1;

}

static int 
check_bounds(PopcountOrdering *popcount_ordering,
	     int *start, int *end, int target_start, int target_end) {
  if (*start > target_end) {
    ordering_no_higher(popcount_ordering);
    //printf("  -- no higher\n");
    return 0;
  }
  if (*end < target_start) {
    ordering_no_lower(popcount_ordering);
    //printf("  -- no lower\n");
    return 0;
  };

  if (*start < target_start) {
    ordering_no_higher(popcount_ordering);
    *start = target_start;
  }
  if (*end > target_end) {
    ordering_no_lower(popcount_ordering);
    *end = target_end;
  }
  return 1;
}

/* count code */
void chemfp_count_tanimoto_arena(
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
    return;
  }
  if (target_start >= target_end) {
    for (query_index = query_start; query_index < query_end; query_index++) {
      /* No possible targets */
      *result_counts++ = 0;
    }
    return;
  }
  if (target_popcount_indicies == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
    query_fp = query_arena + (query_start * query_storage_size);
    for (query_index = query_start; query_index < query_end;
	 query_index++, query_fp += query_storage_size) {
      target_fp = target_arena + (target_start * target_storage_size);
      // Handle the popcount(query) == 0 special case?
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
    return;
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
  } /* go through each of the queries */
}



/**** k-largest code ****/

static int 
klargest_tanimoto_arena_no_popcounts(
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
	int num_allocated,       /* Number of cells allocated */
	int *result_offsets,
	int *result_indicies,
	double *result_scores
				   ) {
  int query_index, target_index;
  int fp_size = (num_bits+7)/8;
  const unsigned char *query_fp, *target_fp;
  double query_threshold, score;
  int heap_size;
  IndexScoreData heap;
  int result_offset=0;


  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {

    if (num_allocated < k) {
      /* Not enough space to store everything, so stop here. */
      return query_index;
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
	  // Since we leave the loop early, I need to advance the pointers
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
    //    printf("Adding %d\n", heap_size);
    //    printf("scores %f %f\n", heap.scores[0], heap.scores[1]);
    result_indicies += heap_size;
    result_scores += heap_size;
    num_allocated -= heap_size;
  } /* Loop through the queries */
  return query_index;
}

int chemfp_klargest_tanimoto_arena(
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
	int num_allocated,       /* Number of cells allocated */
	int *result_offsets,
	int *result_indicies,
	double *result_scores
				   ) {

  int heap_size, fp_size;
  int query_popcount, target_popcount, intersect_popcount;
  double score, best_possible_score, popcount_sum, query_threshold;
  const unsigned char *query_fp, *target_fp;
  int query_index, target_index;
  int start, end;
  PopcountOrdering popcount_ordering;
  IndexScoreData heap;
  int result_offset;

  /* This is C. We don't check for illegal input values. */

  if (query_start >= query_end) {
    return query_start;
  }
  result_offset = 0;
  *result_offsets++ = 0;  // The first offset is always 0

  /* k == 0 is a valid input, and of course the result is no matches */
  if (k == 0) {
    for (query_index = query_start; query_index < query_end; query_index++) {
      *result_offsets++ = 0;
    }
    return query_index;
  }
  fp_size = (num_bits+7)/8;

  if (target_popcount_indicies == NULL) {
    /* precomputed targets aren't available. Use the slower algorithm. */
    return klargest_tanimoto_arena_no_popcounts(
	k, threshold, fp_size,
	query_storage_size, query_arena, query_start, query_end,
	target_storage_size, target_arena, target_start, target_end,
	num_allocated, result_offsets, result_indicies, result_scores);
  }

  /* Loop through the query fingerprints */
  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {
    printf("Query index %d\n", query_index);
    if (num_allocated < k) {
      /* Not enough space to store everything, so stop here. */
      return query_index;
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
    init_popcount_ordering(&popcount_ordering, query_popcount, num_bits);

    heap_size = 0;
    heap.indicies = result_indicies;
    heap.scores = result_scores;

    /* Look through the sections of the arena in optimal popcount order */
    while (next_popcount(&popcount_ordering, query_threshold)) {
      target_popcount = popcount_ordering.popcount;
      best_possible_score = popcount_ordering.score;
      //printf("popcount %d %f\n", target_popcount, best_possible_score);

      /* If we can't beat the query threshold then we're done with the targets */
      if (best_possible_score < query_threshold) {
	break;
      }

      /* Scan through the targets which have the given popcount */
      start = target_popcount_indicies[target_popcount];
      end = target_popcount_indicies[target_popcount+1];
      //printf("  start %d end %d\n", start, end);
      
      if (!check_bounds(&popcount_ordering, &start, &end, target_start, target_end)) {
	continue;
      }
      //printf("  (adjusted) start %d end %d\n", start, end);

      /* Iterate over the target fingerprints */
      target_fp = target_arena + start*target_storage_size;
      popcount_sum = (double)(query_popcount + target_popcount);

      target_index = start;

      /* There are fewer than 'k' elements in the heap*/
      if (heap_size < k) {
	for (; target_index<end; target_index++, target_fp += target_storage_size) {
	  intersect_popcount = chemfp_byte_intersect_popcount(fp_size, query_fp, target_fp);
	  score = intersect_popcount / (popcount_sum - intersect_popcount);
	  //printf("  testing %d score %f from %d %d %d\n", target_index, score,
	  //	 intersect_popcount, query_popcount, target_popcount);

	  /* The heap isn't full; only check if we're at or above the query threshold */
	  if (score >= query_threshold) {
	    heap.indicies[heap_size] = target_index;
	    heap.scores[heap_size] = score;
	    heap_size++;
	    if (heap_size == k) {
	      chemfp_heapq_heapify(heap_size, &heap,  (chemfp_heapq_lt) double_score_lt,
				   (chemfp_heapq_swap) double_score_swap);
	      query_threshold = heap.scores[0];
	      // We're going to jump to the "heap is full" section
	      // Since we leave the loop early, I need to advance the pointers
	      target_index++;
	      target_fp += target_storage_size;
	      //printf("w00t\n");
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

	//printf("  full testing %d score %f from %d %d %d\n", target_index, score,
	//       intersect_popcount, query_popcount, target_popcount);
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
      memcpy(result_indicies, heap.indicies, heap_size * sizeof(int));
      memcpy(result_scores, heap.scores, heap_size * sizeof(double));

      /* Move the pointers so I can report the next results */
      //      printf("Adding %d\n", heap_size);
      //      printf("scores %f %f\n", heap.scores[0], heap.scores[1]);
      result_indicies += heap_size;
      result_scores += heap_size;
      num_allocated -= heap_size;
    }

  } /* looped over all queries */

  return query_index;
}

typedef struct {
  int popcount;
  int index;
} PopcountReorder;

static int compare_by_popcount(const void *left_p, const void *right_p) {
  const PopcountReorder *left = (PopcountReorder *) left_p;
  const PopcountReorder *right = (PopcountReorder *) right_p;
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
	unsigned char *new_arena, int *popcount_indicies) {

  int num_fingerprints, popcount;
  int fp_size = (num_bits+7)/8;
  int i;
  PopcountReorder *reordering;
  const unsigned char *fp;

  if (end <= start) {
    for (i=0; i<=(num_bits+1); i++) {
      *popcount_indicies++ = 0;
    }
    return 0;
  }
  num_fingerprints = end - start;
  reordering = malloc(num_fingerprints * sizeof(PopcountReorder));
  if (!reordering) {
    return CHEMFP_NO_MEM;
  }

  fp = arena + (start * storage_size);
  for (i = start; i<end; i++, fp += storage_size) {
    popcount = chemfp_byte_popcount(fp_size, fp);
    reordering[i].popcount = popcount;
    reordering[i].index = i;
  }
  qsort(reordering, num_fingerprints, sizeof(PopcountReorder), compare_by_popcount);


  /* Build the new arena based on the values in the old arena */

  for (i=0; i<num_fingerprints; i++) {
    memcpy(new_arena, arena+(reordering[i].index * storage_size), storage_size);
    new_arena += storage_size;
  }
  
  /* Create the popcount indicies */
  if (popcount_indicies != NULL) {
    /* Since we've sorted by popcount, this is easy */
    popcount = 0;
    *popcount_indicies++ = 0;
    for (i=0; i<num_fingerprints; i++) {
      while (popcount < reordering[i].popcount) {
	*popcount_indicies++ = i;
	popcount++;
	if (popcount == num_bits) {
	  // We are at or above the limit. We can stop now.
	  i = num_fingerprints+1;
	  break;
	  // Note: with incorrupte data it is possible
	  // that ->popcount can be > num_bits. This is
	  // undefined behavior. I get to do what I want.
	  // I decided to treat them as having "max_popcount" bits.
	  // After all, I don't want corrupt data to crash the
	  // system, and no one is going to validate the input
	  // fingerprints for correctness each time.
	}
      }
    }
    /* Finish up the high end */
    while (popcount <= num_bits) {
      *popcount_indicies++ = i;
      popcount++;
    }
  }

  free(reordering);
  return num_fingerprints;
}
			   
