#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "heapq.h"
#include "chemfp.h"
#include "chemfp_internal.h"

#if defined(_OPENMP)
  #include <omp.h>
#endif

typedef struct {
  int *indices;
  double *scores;
} IndexScoreData;

static int double_score_lt(IndexScoreData *data, int i, int j) {
  if (data->scores[i] < data->scores[j])
    return 1;
  if (data->scores[i] > data->scores[j])
    return 0;
  /* Sort in descending order by index. (XXX important or overkill?) */
  return (data->indices[i] >= data->indices[j]);
}
static void double_score_swap(IndexScoreData *data, int i, int j) {
  int tmp_index = data->indices[i];
  double tmp_score = data->scores[i];
  data->indices[i] = data->indices[j];
  data->scores[i] = data->scores[j];
  data->indices[j] = tmp_index;
  data->scores[j] = tmp_score;
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

        /* Query arena, start and end indices */
        int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size,
        const unsigned char *target_arena, int target_start, int target_end,

        /* Target popcount distribution information */
        int *target_popcount_indices,

        /* Results go into these arrays  */
        int *result_counts
                                   ) {

  int query_index, target_index;
  const unsigned char *query_fp, *target_fp;
  int start, end;
  int count;
  int fp_size = (num_bits+7) / 8;
  double score, popcount_sum;
  int query_popcount, start_target_popcount, end_target_popcount;
  int target_popcount;
  int intersect_popcount;

  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;
  
  if (query_start >= query_end) {
    /* No queries */
    return 0;
  }
  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  /* (Not really needed unless you trap IEE 754 overflow errors) */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    for (query_index = query_start; query_index < query_end; query_index++) {
      /* No possible targets */
      *result_counts++ = 0;
    }
    return query_index-query_start;
  }

  if (threshold <= 0.0) {
    /* Everything will match, so there's no need to figure that out */
    for (query_index = query_start; query_index < query_end; query_index++) {
      *result_counts++ = (target_end-target_start);
    }
    return query_index-query_start;
  }

  if (target_popcount_indices == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
    query_fp = query_arena + (query_start * query_storage_size);
    for (query_index = query_start; query_index < query_end;
         query_index++, query_fp += query_storage_size) {
      target_fp = target_arena + (target_start * target_storage_size);
      /* Handle the popcount(query) == 0 special case? */
      count = 0;

      for (target_index = target_start; target_index < target_end;
           target_index++, target_fp += target_storage_size) {
        score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
        if (score >= threshold) {
          count++;
        }
      }
      *result_counts++ = count;
    }
    return query_index-query_start;
  }
                                                   
  /* Choose popcounts optimized for this case */
  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);

  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the ordering because it's supposed to find everything */

  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {
    
    query_popcount = calc_popcount(fp_size, query_fp);
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
      start_target_popcount = (int)(query_popcount * threshold);
      end_target_popcount = (int)(ceil(query_popcount / threshold));
      if (end_target_popcount > num_bits) {
        end_target_popcount = num_bits;
      }
    }
    count = 0;
    #pragma omp parallel for private(start, end, target_fp, popcount_sum, target_index, intersect_popcount, score) schedule(dynamic)
    for (target_popcount=start_target_popcount; target_popcount<=end_target_popcount;
         target_popcount++) {
      start = target_popcount_indices[target_popcount];
      end = target_popcount_indices[target_popcount+1];
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
        intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
        score = intersect_popcount / (popcount_sum - intersect_popcount);
        if (score >= threshold) {
          #pragma omp atomic
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

        /* Query arena, start and end indices */
        int query_storage_size, const unsigned char *query_arena,
        int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size, const unsigned char *target_arena,
        int target_start, int target_end,

        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go into these arrays  */
        chemfp_threshold_result *results) {

  int query_index, target_index;
  const unsigned char *query_fp, *target_fp;
  int start, end;
  int fp_size = (num_bits+7) / 8;
  double score;
  int query_popcount, start_target_popcount, end_target_popcount;
  int target_popcount;
  int intersect_popcount, popcount_sum;

  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;
  
  if (query_start >= query_end) {
    /* No queries */
    return CHEMFP_OK;
  }

  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  /* (Not really needed unless you trap IEE 754 overflow errors) */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    return CHEMFP_OK;
  }
  if (target_popcount_indices == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
    query_fp = query_arena + (query_start * query_storage_size);
    for (query_index = query_start; query_index < query_end;
         query_index++, query_fp += query_storage_size) {
      target_fp = target_arena + (target_start * target_storage_size);
      /* Handle the popcount(query) == 0 special case? */
      for (target_index = target_start; target_index < target_end;
           target_index++, target_fp += target_storage_size) {
        score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
        if (score >= threshold) {
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, score)) {
            return CHEMFP_NO_MEM;
          };
        }
      }
    }
    return CHEMFP_OK;
  }
  

  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);
  

  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the ordering because it's supposed to find everything */

  query_fp = query_arena + (query_start * query_storage_size);
  for (query_index = query_start; query_index < query_end;
       query_index++, query_fp += query_storage_size) {

    query_popcount = calc_popcount(fp_size, query_fp);
    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
        for (target_index = target_start; target_index < target_end; target_index++) {
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, 0.0)) {
            return CHEMFP_NO_MEM;
          };
        }
      }
      continue;
    }
    /* Figure out which fingerprints to search */
    if (threshold == 0.0) {
      start_target_popcount = 0;
      end_target_popcount = num_bits;
    } else {
      start_target_popcount = (int)(query_popcount * threshold);
      end_target_popcount = (int)(ceil(query_popcount / threshold));
      if (end_target_popcount > num_bits) {
        end_target_popcount = num_bits;
      }
    }

    for (target_popcount=start_target_popcount; target_popcount<=end_target_popcount;
         target_popcount++) {
      start = target_popcount_indices[target_popcount];
      end = target_popcount_indices[target_popcount+1];
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
        intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
        score = ((double) intersect_popcount) / (popcount_sum - intersect_popcount);
        if (score >= threshold) {
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, score)) {
            return CHEMFP_NO_MEM;
          };
        }
      }
    }
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

        /* Query arena, start and end indices */
        int query_storage_size, const unsigned char *query_arena,
        int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size, const unsigned char *target_arena,
        int target_start, int target_end,

        /* Results go into these arrays  */
        int *result_offsets,
        int num_cells,
        int *result_indices,
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
    heap.indices = result_indices;
    heap.scores = result_scores;
    
    target_fp = target_arena + (target_start * query_storage_size);
    target_index = target_start;

    for (; target_index < target_end;
         target_index++, target_fp += target_storage_size) {
      score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
      if (score >= query_threshold) {
        heap.indices[heap_size] = target_index;
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
          heap.indices[0] = target_index;
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
    memcpy(result_indices, heap.indices, heap_size * sizeof(int));
    memcpy(result_scores, heap.scores, heap_size * sizeof(double));

    /* Move the pointers so I can report the next results */
    result_indices += heap_size;
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

        /* Query arena, start and end indices */
        int query_storage_size, const unsigned char *query_arena,
        int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size, const unsigned char *target_arena,
        int target_start, int target_end,

        /* Target popcount distribution information */
        int *target_popcount_indices,

        /* Results go into these arrays  */
        int *result_offsets,
        int num_cells,
        int *result_indices,
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

  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;

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

  if (target_popcount_indices == NULL) {
    /* precomputed targets aren't available. Use the slower algorithm. */
    return knearest_tanimoto_arena_no_popcounts(
        k, threshold, num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        result_offsets, num_cells, result_indices, result_scores);
  }

  /* Choose popcounts optimized for this case */
  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);

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
    query_popcount = calc_popcount(fp_size, query_fp);

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
    heap.indices = result_indices;
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
      start = target_popcount_indices[target_popcount];
      end = target_popcount_indices[target_popcount+1];
      
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
          intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
          score = intersect_popcount / (popcount_sum - intersect_popcount);

          /* The heap isn't full; only check if we're at or above the query threshold */
          if (score >= query_threshold) {
            heap.indices[heap_size] = target_index;
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
        intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
        score = intersect_popcount / (popcount_sum - intersect_popcount);

        /* We need to be strictly *better* than what's in the heap */
        if (score > query_threshold) {
          heap.indices[0] = target_index;
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
      result_indices += heap_size;
      result_scores += heap_size;
      num_cells -= heap_size;
    }

  } /* looped over all queries */

  return query_index-query_start;
}
