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
  }

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


/**** Support for the k-nearest code ****/

static int double_score_lt(chemfp_threshold_result *result, int i, int j) {
  if (result->scores[i] < result->scores[j])
    return 1;
  if (result->scores[i] > result->scores[j])
    return 0;
  /* Sort in descending order by index. (XXX important or overkill?) */
  return (result->indices[i] >= result->indices[j]);
}
static void double_score_swap(chemfp_threshold_result *result, int i, int j) {
  int tmp_index = result->indices[i];
  double tmp_score = result->scores[i];
  result->indices[i] = result->indices[j];
  result->scores[i] = result->scores[j];
  result->indices[j] = tmp_index;
  result->scores[j] = tmp_score;
}


void chemfp_knearest_results_finalize(chemfp_threshold_result *results_start,
                                      chemfp_threshold_result *results_end) {
  chemfp_threshold_result *result;
  for (result = results_start; result < results_end; result++) {
    /* Sort the elements */
    chemfp_heapq_heapsort(result->num_hits, result, (chemfp_heapq_lt) double_score_lt,
                          (chemfp_heapq_swap) double_score_swap);
  }
}
                              
/***** Define the main interface code ***/

#if defined(_OPENMP)

#define RESULT static int
#define RENAME(name) _ ## name ## _single
#define USE_OPENMP 0
#include "search_core.c"
#undef RENAME
#undef USE_OPENMP

#define RENAME(name) _ ## name ## _openmp
#define USE_OPENMP 1
#include "search_core.c"
#undef RENAME
#undef USE_OPENMP
#undef RESULT

/* Dispatch based on the number of threads in use */

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
  if (chemfp_get_num_threads() <= 1)  {
    return _chemfp_count_tanimoto_arena_single(
                           threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, result_counts);
  } else {
    return _chemfp_count_tanimoto_arena_openmp(
                           threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, result_counts);
  }
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

        /* Results go here */
        chemfp_threshold_result *results) {

  if (chemfp_get_num_threads() <= 1) {
    return _chemfp_threshold_tanimoto_arena_single(
                           threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, results);
  } else {
    return _chemfp_threshold_tanimoto_arena_openmp(
                           threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, results);
  }
}

int chemfp_knearest_tanimoto_arena(
        /* Find the 'k' nearest items */
        int k,
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

        /* Results go here */
        chemfp_threshold_result *results) {

  if (chemfp_get_num_threads() <= 1) {
    return _chemfp_knearest_tanimoto_arena_single(
                           k, threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, results);
  } else {
    return _chemfp_knearest_tanimoto_arena_openmp(
                           k, threshold, num_bits,
                           query_storage_size, query_arena, query_start, query_end,
                           target_storage_size, target_arena, target_start, target_end,
                           target_popcount_indices, results);
  }
}


#else

/* Not compiling for OpenMP; don't need the run-time switch */
/* Instead, just rename the function */

#define RESULT int
#define RENAME(name) name
#define USE_OPENMP 0
#include "search_core.c"
#undef USE_OPENMP
#undef RENAME
#undef RESULT

#endif



/***** Special support for the NxN symmetric case ******/

/* TODO: implement the k-nearest variant. It's harder because a k-nearest
   search, combined with the Swamidass and Baldi search limits, is not reflexive. */

#define MAX(x, y) ((x) > (y) ? (x) : (y))

int chemfp_count_tanimoto_hits_arena_symmetric(
        /* Count all matches within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Fingerprint arena */
        int storage_size, const unsigned char *arena,

        /* Row start and end indices */
        int query_start, int query_end,

        /* Column start and end indices */
        int target_start, int target_end,

        /* Target popcount distribution information */
        int *target_popcount_indices,

        /* Results _increment_ existing values in the array - remember to initialize! */
        int *result_counts
                                          ) {
  int fp_size = (num_bits+7) / 8;
  int query_index, target_index;
  int start, end;
  int query_popcount, target_popcount;
  int start_target_popcount, end_target_popcount, intersect_popcount;
  int count;
  double popcount_sum, score;
  const unsigned char *query_fp, *target_fp;
  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;

  /* Check that we're not obviously in the lower triangle */
  if (query_start >= target_end) {  /* No possible hits */
    return CHEMFP_OK;
  }

  /* Shift the target towards the upper triangle, if needed */
  if (target_start < query_start) {
    target_start = query_start;
  }

  /* Check for edge cases */
  if ((query_start >= query_end) ||
      (target_start >= target_end) ||
      (threshold > 1.0)) {
    return CHEMFP_OK;
  }

  if (threshold <= 0.0) {
    /* By definition, everything matches */
    for (query_index = query_start; query_index < query_end; query_index++) {
      start = MAX(query_index+1, target_start);
      end = MAX(query_index+1, target_end);
      if (start < end) {
        result_counts[query_index] += (end - start);
      }
    }
    return CHEMFP_OK;
  }


  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  /* (Not really needed unless you trap IEEE 754 overflow errors) */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }

  /* target_popcount_indices must exist; if you don't care for the factor */
  /* of two performance increase by precomputing/presorting based on popcount */
  /* then why are you interested in the factor of two based on symmetry? */
                                                   
  /* Choose popcount methods optimized for this case */
  calc_popcount = chemfp_select_popcount(num_bits, storage_size, arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, storage_size, arena, storage_size, arena);

  /* This uses the limits from Swamidass and Baldi */
  for (query_index = query_start; query_index < query_end; query_index++) {
    query_fp = arena + (query_index * storage_size);
    query_popcount = calc_popcount(fp_size, query_fp);

    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      continue;
    }
    /* Figure out which fingerprints to search */
    start_target_popcount = (int)(query_popcount * threshold);
    end_target_popcount = (int)(ceil(query_popcount / threshold));
    if (end_target_popcount > num_bits) {
      end_target_popcount = num_bits;
    }

    count = 0;
    for (target_popcount = start_target_popcount; target_popcount <= end_target_popcount;
         target_popcount++) {
      start = target_popcount_indices[target_popcount];
      end = target_popcount_indices[target_popcount+1];
      if (start < target_start) {
        start = target_start;
      }
      start = MAX(query_index+1, start);
      if (end > target_end) {
        end = target_end;
      }

      target_fp = arena + (start * storage_size);
      popcount_sum = query_popcount + target_popcount;
      for (target_index = start; target_index < end;
           target_index++, target_fp += storage_size) {
        intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
        score = intersect_popcount / (popcount_sum - intersect_popcount);
        if (score >= threshold) {
          /* Can accumulate the score for the row */
          count++;
          /* But can't for the symmetric match */
          result_counts[target_index]++;
        }
      }
    }
    result_counts[query_index] += count;
  } /* went through each of the queries */
  return CHEMFP_OK;
}

int chemfp_threshold_tanimoto_arena_symmetric(
        /* Within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Arena */
        int storage_size, const unsigned char *arena,

        /* start and end indices for the rows and columns */
        int query_start, int query_end,
        int target_start, int target_end,
        
        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go here */
        /* NOTE: This must have enough space for all of the fingerprints! */
        chemfp_threshold_result *results) {

  int fp_size = (num_bits+7) / 8;
  int query_index, target_index;
  int start, end;
  const unsigned char *query_fp, *target_fp;
  int query_popcount, target_popcount;
  int start_target_popcount, end_target_popcount;
  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;
  int numerator, denominator, popcount_sum, intersect_popcount;
  double score;
  int add_hit_error = 0;

  /* Check that we're not obviously in the lower triangle */
  if (query_start >= target_end) {  /* No possible hits */
    return CHEMFP_OK;
  }

  /* Shift the target towards the upper triangle, if needed */
  if (target_start < query_start) {
    target_start = query_start;
  }

  /* Corner cases where I don't need to do anything */
  if ((query_start >= query_end) ||
      (target_start >= target_end) ||
      (threshold < 0)) {
    return CHEMFP_OK;
  }

  /* if (threshold == 0.0) { */ /* TODO: Optimize this case */


  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }
  if (threshold > 1.0) {
    return CHEMFP_OK;
  }

  calc_popcount = chemfp_select_popcount(num_bits, storage_size, arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, storage_size, arena, storage_size, arena);
  
  denominator = num_bits * 10;
  numerator = (int)(threshold * denominator);

  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the search ordering because it's supposed to find everything */
  
  for (query_index = query_start; query_index < query_end; query_index++) {
    query_fp = arena + (query_index * storage_size);
    query_popcount = calc_popcount(fp_size, query_fp);

    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
        /* Only populate the upper triangle */
        target_index = MAX(query_index+1, target_start);
        for (;target_index < target_end; target_index++) {
          if (!_chemfp_add_hit(results+query_index, target_index, 0.0)) {
            add_hit_error = 1;
          }
        }
      }
      continue;
    }
    /* Figure out which fingerprints to search, based on the popcount */
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

      popcount_sum = query_popcount + target_popcount;
      for (target_index = MAX(query_index+1, start); target_index < end; target_index++) {
        target_fp = arena + (target_index * storage_size);
        intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);

        score = ((double) intersect_popcount) / (popcount_sum - intersect_popcount);
        if (denominator * intersect_popcount  >=
            numerator * (popcount_sum - intersect_popcount)) {
          /* Add to the upper triangle */
          if (!_chemfp_add_hit(results+query_index, target_index, score)) {
            add_hit_error = 1;
          }
        }
      }
    }
  } /* went through each of the queries */
  if (add_hit_error) {
    return CHEMFP_NO_MEM;
  }
  return CHEMFP_OK;
}

