#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "chemfp.h"
#include "chemfp_internal.h"

/****************************** Start of TimSort code ************************/

/* The following was modified to support sorting two parallel arrays.
   I also removed the macros. */

/* License

All code [in the repository at https://github.com/swenson/sort],
unless otherwise specified, is hereby licensed under the MIT Public
License:

  Copyright (c) 2010 Christopher Swenson

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Modified by Andrew Dalke. I extracted the TimSort code and modified it
to work with two parallel arrays.

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>


#define SORT_TYPE1 double
#define SORT_TYPE2 int

typedef int (*hit_compare_func)(SORT_TYPE1 left_score, SORT_TYPE1 right_score,
                                SORT_TYPE2 left_index, SORT_TYPE2 right_index);
#define SORT_CMP(x1, y1, x2, y2) hit_compare(x1, y1, x2, y2)

#ifndef CLZ
#ifdef __GNUC__
#define CLZ __builtin_clzll
#else

/* adapted from Hacker's Delight */
static int clzll(uint64_t x) {
  int n;

  if (x == 0) return(64);
  n = 0;
  if (x <= 0x00000000FFFFFFFFL) {n = n + 32; x = x << 32;}
  if (x <= 0x0000FFFFFFFFFFFFL) {n = n + 16; x = x << 16;}
  if (x <= 0x00FFFFFFFFFFFFFFL) {n = n + 8; x = x << 8;}
  if (x <= 0x0FFFFFFFFFFFFFFFL) {n = n + 4; x = x << 4;}
  if (x <= 0x3FFFFFFFFFFFFFFFL) {n = n + 2; x = x << 2;}
  if (x <= 0x7FFFFFFFFFFFFFFFL) {n = n + 1;}
  return n;
}

#define CLZ clzll
#endif
#endif


#define SORT_SWAP1(x,y) ({SORT_TYPE1 __SORT_SWAP_t = (x); (x) = (y); (y) = __SORT_SWAP_t;})
#define SORT_SWAP2(x,y) ({SORT_TYPE2 __SORT_SWAP_t = (x); (x) = (y); (y) = __SORT_SWAP_t;})

#ifndef MAX
#define MAX(x,y) (((x) > (y) ? (x) : (y)))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y) ? (x) : (y)))
#endif

typedef struct {
  int64_t start;
  int64_t length;  
} hits_tim_sort_run_t;


void hits_binary_insertion_sort(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t size, hit_compare_func hit_compare);
void hits_tim_sort(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t size, hit_compare_func hit_compare);


/* Function used to do a binary search for binary insertion sort */
static inline int64_t binary_insertion_find(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const SORT_TYPE1 x1, const SORT_TYPE2 x2, const size_t size, hit_compare_func hit_compare)
{
  int64_t l, c, r;
  l = 0;
  r = size - 1;
  c = r >> 1;
  SORT_TYPE1 lx1, cx1, rx1;
  SORT_TYPE2 lx2, cx2, rx2;
  lx1 = dst1[l];
  lx2 = dst2[l];
  
  /* check for beginning conditions */
  if (SORT_CMP(x1, lx1, x2, lx2) < 0)
    return 0;
  else if (SORT_CMP(x1, lx1, x2, lx2) == 0)
  {
    int64_t i = 1;
    while (SORT_CMP(x1, dst1[i], x2, dst2[i]) == 0) i++;
    return i;
  }
  
  rx1 = dst1[r];
  rx2 = dst2[r];
  /* guaranteed not to be >= rx */
  cx1 = dst1[c];
  cx2 = dst2[c];
  while (1)
  {
    const int val = SORT_CMP(x1, cx1, x2, cx2);
    if (val < 0)
    {
      if (c - l <= 1) return c;
      r = c;
      rx1 = cx1;
      rx2 = cx2;
    }
    else if (val > 0)
    {
      if (r - c <= 1) return c + 1;
      l = c;
      lx1 = cx2;
      lx2 = cx2;
    }
    else
    {
      do
      {
        cx1 = dst1[c];
        cx2 = dst2[++c];
      } while (SORT_CMP(x1, cx1, x2, cx2) == 0);
      return c;
    }
    c = l + ((r - l) >> 1);
    cx1 = dst1[c];
    cx2 = dst2[c];
  }
}

/* Binary insertion sort, but knowing that the first "start" entries are sorted.  Used in timsort. */
static inline void binary_insertion_sort_start(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t start, const size_t size, hit_compare_func hit_compare)
{
  int64_t i;
  for (i = start; i < size; i++)
  {
    int64_t j;
    /* If this entry is already correct, just move along */
    if (SORT_CMP(dst1[i - 1], dst1[i], dst2[i - 1], dst2[i]) <= 0) continue;
    
    /* Else we need to find the right place, shift everything over, and squeeze in */
    SORT_TYPE1 x1 = dst1[i];
    SORT_TYPE2 x2 = dst2[i];
    int64_t location = binary_insertion_find(dst1, dst2, x1, x2, i, hit_compare);
    for (j = i - 1; j >= location; j--)
    {
      dst1[j + 1] = dst1[j];
      dst2[j + 1] = dst2[j];
    }
    dst1[location] = x1;
    dst2[location] = x2;
  }
}

/* Binary insertion sort */
void hits_binary_insertion_sort(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t size, hit_compare_func hit_compare)
{
  binary_insertion_sort_start(dst1, dst2, 1, size, hit_compare);
}

/* timsort implementation, based on timsort.txt */

static inline void reverse_elements(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, int64_t start, int64_t end)
{
  while (1)
  {
    if (start >= end) return;
    SORT_SWAP1(dst1[start], dst1[end]);
    SORT_SWAP2(dst2[start], dst2[end]);
    start++;
    end--;
  }
}

static inline int64_t count_run(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const int64_t start, const size_t size, hit_compare_func hit_compare)
{
  if (size - start == 1) return 1;
  if (start >= size - 2)
  {
    if (SORT_CMP(dst1[size - 2], dst1[size - 1], dst2[size - 2], dst2[size - 1]) > 0)
    {
      SORT_SWAP1(dst1[size - 2], dst2[size - 1]);
      SORT_SWAP2(dst2[size - 2], dst2[size - 1]);
    }
    return 2;
  }
  
  int64_t curr = start + 2;
  
  if (SORT_CMP(dst1[start], dst1[start + 1], dst2[start], dst2[start + 1]) <= 0)
  {
    /* increasing run */
    while (1)
    {
      if (curr == size - 1) break;
      if (SORT_CMP(dst1[curr - 1], dst1[curr], dst2[curr - 1], dst2[curr]) > 0) break;
      curr++;
    }
    return curr - start;    
  }
  else
  {
    /* decreasing run */
    while (1)
    {
      if (curr == size - 1) break;
      if (SORT_CMP(dst1[curr - 1], dst1[curr], dst2[curr - 1], dst2[curr]) <= 0) break;
      curr++;
    }
    /* reverse in-place */
    reverse_elements(dst1, dst2, start, curr - 1);
    return curr - start;
  }
}

static inline int compute_minrun(const uint64_t size)
{
  const int top_bit = 64 - CLZ(size);
  const int shift = MAX(top_bit, 6) - 6;
  const int minrun = size >> shift;
  const uint64_t mask = (1ULL << shift) - 1;
  if (mask & size) return minrun + 1;
  return minrun;
}

#define PUSH_NEXT() do {\
len = count_run(dst1, dst2, curr, size, hit_compare);   \
run = minrun;\
if (run < minrun) run = minrun;\
if (run > size - curr) run = size - curr;\
if (run > len)\
{\
  binary_insertion_sort_start(&dst1[curr], &dst2[curr], len, run, hit_compare);     \
  len = run;\
}\
run_stack[stack_curr++] = (hits_tim_sort_run_t) {curr, len};\
curr += len;\
if (curr == size)\
{\
  /* finish up */ \
  while (stack_curr > 1) \
  { \
    tim_sort_merge(dst1, dst2, run_stack, stack_curr, store, hit_compare);          \
    run_stack[stack_curr - 2].length += run_stack[stack_curr - 1].length; \
    stack_curr--; \
  } \
  if (store->storage1 != NULL)\
  {\
    free(store->storage1);\
    store->storage1 = NULL;\
    free(store->storage2);\
    store->storage2 = NULL;\
  }\
  return;\
}\
}\
while (0)
  
static inline int check_invariant(hits_tim_sort_run_t *stack, const int stack_curr)
{
  if (stack_curr < 2) return 1;
  if (stack_curr == 2)
  {
    const int64_t A = stack[stack_curr - 2].length;
    const int64_t B = stack[stack_curr - 1].length;
    if (A <= B) return 0;
    return 1;
  }
  const int64_t A = stack[stack_curr - 3].length;
  const int64_t B = stack[stack_curr - 2].length;
  const int64_t C = stack[stack_curr - 1].length;
  if ((A <= B + C) || (B <= C)) return 0;
  return 1;
}

typedef struct {
  size_t alloc;
  SORT_TYPE1 *storage1;
  SORT_TYPE2 *storage2;
} hits_temp_storage_t;
  

static inline void tim_sort_resize(hits_temp_storage_t *store, const size_t new_size)
{
  if (store->alloc < new_size)
  {
    SORT_TYPE1 *tempstore1 = realloc(store->storage1, new_size * sizeof(SORT_TYPE1));
    SORT_TYPE2 *tempstore2 = realloc(store->storage2, new_size * sizeof(SORT_TYPE2));
    if (tempstore1 == NULL)
    {
      fprintf(stderr, "Error allocating temporary storage for tim sort: need %lu bytes", sizeof(SORT_TYPE1) * new_size);
      exit(1);
    }
    if (tempstore2 == NULL)
    {
      fprintf(stderr, "Error allocating temporary storage for tim sort: need %lu bytes", sizeof(SORT_TYPE2) * new_size);
      exit(1);
    }
    store->storage1 = tempstore1;
    store->storage2 = tempstore2;
    store->alloc = new_size;
  }
}

static inline void tim_sort_merge(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const hits_tim_sort_run_t *stack, const int stack_curr, hits_temp_storage_t *store, hit_compare_func hit_compare)
{
  const int64_t A = stack[stack_curr - 2].length;
  const int64_t B = stack[stack_curr - 1].length;
  const int64_t curr = stack[stack_curr - 2].start;

  tim_sort_resize(store, MIN(A, B));
  SORT_TYPE1 *storage1 = store->storage1;
  SORT_TYPE2 *storage2 = store->storage2;
  
  int64_t i, j, k;

  /* left merge */
  if (A < B)
  {
    memcpy(storage1, &dst1[curr], A * sizeof(SORT_TYPE1));
    memcpy(storage2, &dst2[curr], A * sizeof(SORT_TYPE2));
    i = 0;
    j = curr + A;
    
    for (k = curr; k < curr + A + B; k++)
    {
      if ((i < A) && (j < curr + A + B))
      {
        if (SORT_CMP(storage1[i], dst1[j], storage2[i], dst2[j]) <= 0) {
          dst1[k] = storage1[i];
          dst2[k] = storage2[i++];
        } else { 
          dst1[k] = dst1[j];
          dst2[k] = dst2[j++];
        }
      }
      else if (i < A)
      {
        dst1[k] = storage1[i];
        dst2[k] = storage2[i++];
      }
      else
      { 
        dst1[k] = dst1[j];
        dst2[k] = dst2[j++];
      }
    }
  }
  /* right merge */
  else
  {    
    memcpy(storage1, &dst1[curr + A], B * sizeof(SORT_TYPE1));
    memcpy(storage2, &dst2[curr + A], B * sizeof(SORT_TYPE2));
    i = B - 1;
    j = curr + A - 1;
    
    for (k = curr + A + B - 1; k >= curr; k--)
    {
      if ((i >= 0) && (j >= curr))
      {
          if (SORT_CMP(dst1[j], storage1[i], dst2[j], storage2[i]) > 0)
          {
            dst1[k] = dst1[j];
            dst2[k] = dst2[j--];
          }
          else
          {
            dst1[k] = storage2[i];
            dst2[k] = storage2[i--];
          }
      }
      else if (i >= 0)
      {
        dst1[k] = storage1[i];
        dst2[k] = storage2[i--];
      }
      else
      {
        dst1[k] = dst1[j];
        dst2[k] = dst2[j--];
      }
    }
  }
}

static inline int tim_sort_collapse(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, hits_tim_sort_run_t *stack, int stack_curr, hits_temp_storage_t *store, const size_t size, hit_compare_func hit_compare)
{
  while (1)
  {
    /* if the stack only has one thing on it, we are done with the collapse */
    if (stack_curr <= 1) break;
    /* if this is the last merge, just do it */
    if ((stack_curr == 2) && (stack[0].length + stack[1].length == size))
    {
      tim_sort_merge(dst1, dst2, stack, stack_curr, store, hit_compare);
      stack[0].length += stack[1].length;
      stack_curr--;
      break;
    }
    /* check if the invariant is off for a stack of 2 elements */
    else if ((stack_curr == 2) && (stack[0].length <= stack[1].length))
    {
      tim_sort_merge(dst1, dst2, stack, stack_curr, store, hit_compare);
      stack[0].length += stack[1].length;
      stack_curr--;
      break;
    }
    else if (stack_curr == 2)
      break;
      
    const int64_t A = stack[stack_curr - 3].length;
    const int64_t B = stack[stack_curr - 2].length;
    const int64_t C = stack[stack_curr - 1].length;
    
    /* check first invariant */
    if (A <= B + C)
    {
      if (A < C)
      {
        tim_sort_merge(dst1, dst2, stack, stack_curr - 1, store, hit_compare);
        stack[stack_curr - 3].length += stack[stack_curr - 2].length;
        stack[stack_curr - 2] = stack[stack_curr - 1];
        stack_curr--;
      }
      else
      {
        tim_sort_merge(dst1, dst2, stack, stack_curr, store, hit_compare);
        stack[stack_curr - 2].length += stack[stack_curr - 1].length;
        stack_curr--;
      }
    }
    /* check second invariant */
    else if (B <= C)
    {
      tim_sort_merge(dst1, dst2, stack, stack_curr, store, hit_compare);
      stack[stack_curr - 2].length += stack[stack_curr - 1].length;
      stack_curr--;      
    }
    else
      break;
  }
  return stack_curr;
}

void hits_tim_sort(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t size, hit_compare_func hit_compare)
{
  if (size < 64)
  {
    hits_binary_insertion_sort(dst1, dst2, size, hit_compare);
    return;
  }
  
  /* compute the minimum run length */
  const int minrun = compute_minrun(size);
  
  /* temporary storage for merges */
  hits_temp_storage_t _store, *store = &_store;
  store->alloc = 0;
  store->storage1 = NULL;
  store->storage2 = NULL;
  
  hits_tim_sort_run_t run_stack[128];
  int stack_curr = 0;
  int64_t len, run;
  int64_t curr = 0;
  
  PUSH_NEXT();
  PUSH_NEXT();
  PUSH_NEXT();
  
  while (1)
  {
    if (!check_invariant(run_stack, stack_curr))
    {
      stack_curr = tim_sort_collapse(dst1, dst2, run_stack, stack_curr, store, size, hit_compare);
      continue;
    }
    PUSH_NEXT();
  }
}



/* heap sort: based on wikipedia */

static inline void heap_sift_down(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const int64_t start, const int64_t end, hit_compare_func hit_compare)
{
  int64_t root = start;
  
  while ((root << 1) <= end)
  {
    int64_t child = root << 1;
    if ((child < end) && (SORT_CMP(dst1[child], dst1[child + 1], dst2[child], dst2[child + 1]) < 0))
      child++;
    if (SORT_CMP(dst1[root], dst1[child], dst2[root], dst2[child]) < 0)
    {
      SORT_SWAP1(dst1[root], dst1[child]);
      SORT_SWAP2(dst2[root], dst2[child]);
      root = child;
    }
    else
      return;
  }
}

static inline void heapify(SORT_TYPE1 *dst1, SORT_TYPE2 *dst2, const size_t size, hit_compare_func hit_compare)
{
  int64_t start = size >> 1;
  while (start >= 0)
  {
    heap_sift_down(dst1, dst2, start, size - 1, hit_compare);
    start--;
  }
}

/****************************** End of TimSort code ************************/

chemfp_search_result *chemfp_alloc_search_results(int size) {
  /* Initializes all of the counts to 0 (and makes everything else nice too) */
  return (chemfp_search_result *) calloc(size, sizeof(chemfp_search_result));
}

void chemfp_free_results(int num_results, chemfp_search_result *results) {
  int i;
  for (i=0; i<num_results; i++) {
    if (results[i].num_hits) {
      free(results[i].scores);
    }
  }
  free(results);
}

int chemfp_get_num_hits(chemfp_search_result *result) {
  return result->num_hits;
}

int chemfp_add_hit(chemfp_search_result *result,
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

int chemfp_search_result_get_hits(chemfp_search_result *result,
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

int chemfp_fill_lower_triangle(int n, chemfp_search_result *results) {
  int i, j;
  int *sizes = (int *) malloc(n * sizeof(int));
  int retval;

  if (!sizes) {
    return CHEMFP_NO_MEM;
  }
  /* Save all of the count information */
  for (i=0; i<n; i++) {
    sizes[i] = chemfp_get_num_hits(results+i);
  }

  retval = CHEMFP_OK;
  for (i=0; i<n; i++) {
    for (j=0; j<sizes[i]; j++) {
      if (!chemfp_add_hit(results+results[i].indices[j], i, results[i].scores[j])) {
        retval = CHEMFP_NO_MEM;
        goto done;
      }
    }
  }

 done:
  free(sizes);
  return retval;
}

typedef struct {
  const char *name;
  hit_compare_func hit_compare;
} sort_method_t;


static int compare_decreasing(double score1, double score2, int index1, int index2) {
  printf("Decreasing %f %f %d %d\n", score1, score2, index1, index2);
  if (score1 < score2) {
    return 1;
  }
  if (score1 > score2) {
    return -1;
  }
  if (index1 < index2) {
    return -1;
  }
  if (index1 > index2) {
    return 1;
  }
  return 0;
}

static int compare_increasing(double score1, double score2, int index1, int index2) {
  printf("Increasing %f %f %d %d\n", score1, score2, index1, index2);
  if (score1 < score2) {
    return -1;
  }
  if (score1 > score2) {
    return 1;
  }
  if (index1 < index2) {
    return -1;
  }
  if (index1 > index2) {
    return 1;
  }
  return 0;
}

sort_method_t sort_methods[] = {
  {"increasing", compare_increasing},
  {"decreasing", compare_decreasing},
  //  {"closest-first", closest_first},
  {NULL}
};

static sort_method_t *chemfp_get_sort_method(const char *name) {
  int i=0;
  while (sort_methods[i].name != NULL) {
    if (!strcmp(name, sort_methods[i].name)) {
      return &sort_methods[i];
    }
    i++;
  }
  return NULL;
};

int chemfp_search_results_sort(int num_results, chemfp_search_result *results,
                               const char *order) {
  int i, num_hits;
  sort_method_t *sort_method = chemfp_get_sort_method(order);
  if (sort_method == NULL) {
    return CHEMFP_BAD_ARG;
  }
  for (i=0; i<num_results; i++) {
    num_hits = results[i].num_hits;
    if (num_hits > 1) {
      hits_tim_sort(results[i].scores, results[i].indices, num_hits,
                    sort_method->hit_compare);
    }
  }
  return CHEMFP_OK;
}

int chemfp_search_result_sort(chemfp_search_result *result, const char *order) {
  int num_hits;
  sort_method_t *sort_method = chemfp_get_sort_method(order);
  if (sort_method == NULL) {
    return CHEMFP_BAD_ARG;
  }
  num_hits = result->num_hits;
  if (num_hits > 1) {
      hits_tim_sort(result->scores, result->indices, num_hits,
                    sort_method->hit_compare);
  }
  return CHEMFP_OK;
}

void chemfp_search_result_clear(chemfp_search_result *result) {
  if (result->num_hits != 0) {
    result->num_hits=0;
    free(result->scores);
    result->scores = NULL;
    result->indices = NULL;
  }
}
