/* This is a rather cumbersome solution to two problems I have with OpenMP.

 1) multiple threads and OpenMP don't mix on a Mac. It segfaults
 during the first openmp call. I want people to be able to use chemfp
 in multi-threaded environments, even with diminished performance, so
 the single thread version should not go through the OpenMP path.

 2) I measured a roughly 5% performance penalty hit with a single
 thread using OpenMP vs. the code compiled without OpenMP.

My solution is to compile the core code twice, one for each path. The
RENAME macro rewrites

   RESULT RENAME(chemfp_count_tanimoto_arena)

to one of:

 static int _chemfp_count_tanimoto_arena_single -- single-threaded, compiler supports OpenMP
 static int _chemfp_count_tanimoto_arena_openmp -- multiple OpenMP threads
 int chemfp_count_tanimoto_arena -- single-threaded, compiler does not support OpenMP

depending on the circumstances. In a normal build, where OpenMP is
available, then this file will be #include'd twice.

*/

/* count code */
RESULT RENAME(chemfp_count_tanimoto_arena)(
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
    return CHEMFP_OK;
  }
  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  /* (Not really needed unless you trap IEEE 754 overflow errors) */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    for (query_index = 0; query_index < (query_end-query_start); query_index++) {
      /* No possible targets */
      result_counts[query_index] = 0;
    }
    return CHEMFP_OK;
  }

  if (threshold <= 0.0) {
    /* Everything will match, so there's no need to figure that out */
    for (query_index = 0; query_index < (query_end-query_start); query_index++) {
      result_counts[query_index] = (target_end - target_start);
    }
    return CHEMFP_OK;
  }

  if (target_popcount_indices == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
#if USE_OPENMP == 1
    #pragma omp parallel for private(query_fp, target_fp, count, target_index, score) schedule(dynamic)
#endif
    for (query_index = 0; query_index < (query_end-query_start); query_index++) {
      query_fp = query_arena + (query_start + query_index) * query_storage_size;
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
      result_counts[query_index] = count;
    }
    return CHEMFP_OK;
  }
                                                   
  /* Choose popcounts optimized for this case */
  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);

  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the search ordering because it's supposed to find everything */
#if USE_OPENMP == 1
  #pragma omp parallel for \
      private(query_fp, query_popcount, start_target_popcount, end_target_popcount, \
          count, target_popcount, start, end, target_fp, popcount_sum, target_index, intersect_popcount, score) \
      schedule(dynamic)
#endif
  for (query_index = 0; query_index < (query_end-query_start); query_index++) {
    query_fp = query_arena + (query_start + query_index) * query_storage_size;
    query_popcount = calc_popcount(fp_size, query_fp);
    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
        result_counts[query_index] = (target_end - target_start);
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
    for (target_popcount = start_target_popcount; target_popcount <= end_target_popcount;
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
          count++;
        }
      }
    }
    result_counts[query_index] = count;
  } /* went through each of the queries */
  return CHEMFP_OK;
}

RESULT RENAME(chemfp_threshold_tanimoto_arena)(
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

  int query_index, target_index;
  const unsigned char *query_fp, *target_fp;
  int start, end;
  int fp_size = (num_bits+7) / 8;
  double score;
  int query_popcount, start_target_popcount, end_target_popcount;
  int target_popcount;
  int intersect_popcount, popcount_sum;
  int numerator, denominator;
  int add_hit_error = 0;

  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;
  
  if (query_start >= query_end) {
    /* No queries */
    return CHEMFP_OK;
  }

  /* Prevent overflow if someone uses a threshold of, say, 1E-80 */
  /* (Not really needed unless you trap IEEE 754 overflow errors) */
  if (threshold > 0.0 && threshold < 1.0/num_bits) {
    threshold = 0.5 / num_bits;
  }
  if ((target_start >= target_end) || threshold > 1.0) {
    return CHEMFP_OK;
  }

  if (target_popcount_indices == NULL) {
    /* Handle the case when precomputed targets aren't available. */
    /* This is a slower algorithm because it tests everything. */
#if USE_OPENMP == 1
    #pragma omp parallel for private(query_fp, target_fp, target_index, score) schedule(dynamic)
#endif
    for (query_index = query_start; query_index < query_end; query_index++) {
      query_fp = query_arena + (query_index * query_storage_size);
      target_fp = target_arena + (target_start * target_storage_size);
      /* Handle the popcount(query) == 0 special case? */
      for (target_index = target_start; target_index < target_end;
           target_index++, target_fp += target_storage_size) {
        score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
        if (score >= threshold) {
#if USE_OPENMP == 1
          #pragma omp critical (add_hit_threshold)
#endif
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, score)) {
            add_hit_error = 1;
          }
        }
      }
    }
    if (add_hit_error) {
      return CHEMFP_NO_MEM;
    }
    return CHEMFP_OK;
  }
  

  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);
  
  denominator = num_bits * 10;
  numerator = (int)(threshold * denominator);

  /* This uses the limits from Swamidass and Baldi */
  /* It doesn't use the search ordering because it's supposed to find everything */

#if USE_OPENMP == 1
  #pragma omp parallel for \
      private(query_fp, query_popcount, target_index, target_fp, start_target_popcount, \
          end_target_popcount, target_popcount, start, end, popcount_sum, intersect_popcount, score) \
      schedule(dynamic)
#endif
  for (query_index = query_start; query_index < query_end; query_index++) {
    query_fp = query_arena + (query_index * query_storage_size);
    query_popcount = calc_popcount(fp_size, query_fp);

    /* Special case when popcount(query) == 0; everything has a score of 0.0 */
    if (query_popcount == 0) {
      if (threshold == 0.0) {
#if USE_OPENMP == 1
        #pragma omp critical (add_hit_threshold)
#endif
        for (target_index = target_start; target_index < target_end; target_index++) {
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, 0.0)) {
            add_hit_error = 1;
          }
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

        /* In my timings (on a Mac), the comparison against a double was a hotspot, */
        /* but division is not. I switch to integer math and gained a 3-4% performance, */
        /* at the cost of slightly more complicated code. */
        score = ((double) intersect_popcount) / (popcount_sum - intersect_popcount);
        if (denominator * intersect_popcount  >=
            numerator * (popcount_sum - intersect_popcount)) {
#if USE_OPENMP == 1
          #pragma omp critical (add_hit_threshold)
#endif
          if (!_chemfp_add_hit(results+(query_index-query_start), target_index, score)) {
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



static int 
RENAME(knearest_tanimoto_arena_no_popcounts)(
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
        chemfp_threshold_result *results
                                   ) {
  int query_index, target_index;
  int fp_size = (num_bits+7)/8;
  const unsigned char *query_fp, *target_fp;
  double query_threshold, score;
  chemfp_threshold_result *result;

  for (query_index = 0; query_index < (query_end-query_start); query_index++) {
    query_fp = query_arena + (query_start+query_index) * query_storage_size;

    result = results+query_index;
    query_threshold = threshold;
    
    target_fp = target_arena + (target_start * query_storage_size);
    target_index = target_start;

    for (; target_index < target_end;
         target_index++, target_fp += target_storage_size) {
      score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);
      if (score >= query_threshold) {
        _chemfp_add_hit(result, target_index, score);
        if (result->num_hits == k) {
          chemfp_heapq_heapify(k, result, (chemfp_heapq_lt) double_score_lt,
                               (chemfp_heapq_swap) double_score_swap);
          query_threshold = result->scores[0];
          /* Since we leave the loop early, I need to advance the pointers */
          target_index++;
          target_fp += target_storage_size;
          break;
        }
      }
    }
    /* Either we've reached the end of the fingerprints or the heap is full */
    if (result->num_hits == k) {
      /* Continue scanning through the fingerprints */
      for (; target_index < target_end;
           target_index++, target_fp += target_storage_size) {
        score = chemfp_byte_tanimoto(fp_size, query_fp, target_fp);

        /* We need to be strictly *better* than what's in the heap */
        if (score > query_threshold) {
          result->indices[0] = target_index;
          result->scores[0] = score;
          chemfp_heapq_siftup(k, result, 0, (chemfp_heapq_lt) double_score_lt,
                              (chemfp_heapq_swap) double_score_swap);
          query_threshold = result->scores[0];
        } /* heapreplaced the old smallest item with the new item */
      }
      /* End of the fingerprint scan */
    } else {
      /* The heap isn't full, so we haven't yet heapified it. */
      chemfp_heapq_heapify(result->num_hits, result,  (chemfp_heapq_lt) double_score_lt,
                           (chemfp_heapq_swap) double_score_swap);
    }
  } /* Loop through the queries */

  return query_index-query_start;
}


RESULT RENAME(chemfp_knearest_tanimoto_arena)(
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
        chemfp_threshold_result *results
                                   ) {

  int fp_size;
  int query_popcount, target_popcount, intersect_popcount;
  double score, best_possible_score, popcount_sum, query_threshold;
  const unsigned char *query_fp, *target_fp;
  int query_index, target_index;
  int start, end;
  PopcountSearchOrder popcount_order;
  chemfp_threshold_result *result;
  
  chemfp_popcount_f calc_popcount;
  chemfp_intersect_popcount_f calc_intersect_popcount;

  /* This is C. We don't check for illegal input values. */

  if (query_start >= query_end) {
    return 0;
  }
  /* k == 0 is a valid input, and of course the result is no matches */
  if (k == 0) {
    return CHEMFP_OK;
  }
  fp_size = (num_bits+7)/8;

  if (target_popcount_indices == NULL) {
    /* precomputed targets aren't available. Use the slower algorithm. */
    return RENAME(knearest_tanimoto_arena_no_popcounts)(
        k, threshold, num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        results);
  }

  /* Choose popcounts optimized for this case */
  calc_popcount = chemfp_select_popcount(num_bits, query_storage_size, query_arena);
  calc_intersect_popcount = chemfp_select_intersect_popcount(
                num_bits, query_storage_size, query_arena,
                target_storage_size, target_arena);

  /* Loop through the query fingerprints */
  for (query_index=0; query_index < (query_end-query_start); query_index++) {
    result = results+query_index;
    query_fp = query_arena + (query_start+query_index) * query_storage_size;

    query_threshold = threshold;
    query_popcount = calc_popcount(fp_size, query_fp);

    if (query_popcount == 0) {
      /* By definition this will never return hits. Even if threshold == 0.0. */
      /* (I considered returning the first k hits, but that's chemically meaningless.) */
      /* XXX change this. Make it returns the first k hits */
      continue;
    }

    /* Search the bins using the ordering from Swamidass and Baldi.*/
    init_search_order(&popcount_order, query_popcount, num_bits);

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
      if (result->num_hits < k) {
        for (; target_index<end; target_index++, target_fp += target_storage_size) {
          intersect_popcount = calc_intersect_popcount(fp_size, query_fp, target_fp);
          score = intersect_popcount / (popcount_sum - intersect_popcount);

          /* The heap isn't full; only check if we're at or above the query threshold */
          if (score >= query_threshold) {
            _chemfp_add_hit(result, target_index, score);
            if (result->num_hits == k) {
              chemfp_heapq_heapify(k, result,  (chemfp_heapq_lt) double_score_lt,
                                   (chemfp_heapq_swap) double_score_swap);
              query_threshold = result->scores[0];
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
          result->indices[0] = target_index;
          result->scores[0] = score;
          chemfp_heapq_siftup(k, result, 0, (chemfp_heapq_lt) double_score_lt,
                              (chemfp_heapq_swap) double_score_swap);
          query_threshold = result->scores[0];
          if (query_threshold >= best_possible_score) {
            /* we can't do any better in this section (or in later ones) */
            break;
          }
        } /* heapreplaced the old smallest item with the new item */
      } /* looped over fingerprints */
    } /* Went through all the popcount regions */

    /* We have scanned all the fingerprints. Is the heap full? */
    if (result->num_hits < k) {
      /* Not full, so need to heapify it. */
      chemfp_heapq_heapify(result->num_hits, result, (chemfp_heapq_lt) double_score_lt,
                           (chemfp_heapq_swap) double_score_swap);
    }
  } /* looped over all queries */
  return CHEMFP_OK;
}

