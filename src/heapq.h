/**** Low-level heap operations, for the best-of-N algorithms ****/

/* These are internal data types and functions. While they may */
/* be available in the library, do not call them directly. */


/* Compare two items in the heap. Return -1 on error, 1 for lt, otherwise 0 */
typedef int (*chemfp_heapq_lt)(void *data, int i, int j);

/* Swap two items in the heap. This function must never fail. */
typedef void (*chemfp_heapq_swap)(void *data, int i, int j);

/* Call after replacing the first element in a heapified list */
int chemfp_heapq_siftup(int len, void *heap, int pos,
                        chemfp_heapq_lt lt, chemfp_heapq_swap swap);

/* Convert the un-ordered list into a heap */
int chemfp_heapq_heapify(int len, void *heap,
                         chemfp_heapq_lt lt, chemfp_heapq_swap swap);

/* Must heapify first */
int chemfp_heapq_heapsort(int len, void *heap,
                          chemfp_heapq_lt lt, chemfp_heapq_swap swap);
