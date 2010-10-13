/* Returns -1 on error, 1 for lt, otherwise 0 */
typedef int (*chemfp_heapq_lt)(void *data, int i, int j);

/* Must never fail */
typedef void (*chemfp_heapq_swap)(void *data, int i, int j);

int chemfp_heapq_heapify(int len, void *heap,
						 chemfp_heapq_lt lt, chemfp_heapq_swap swap);

/* Call after replacing the first element in a heapified list */
int chemfp_heapq_siftup(int len, void *heap, int pos,
						chemfp_heapq_lt lt, chemfp_heapq_swap swap);

/* Must heapify first */
int chemfp_heapq_heapsort(int len, void *heap,
						  chemfp_heapq_lt lt, chemfp_heapq_swap swap);
