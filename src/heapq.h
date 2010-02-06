/* Returns -1 on error, 1 for lt, otherwise 0 */
typedef int (*heapq_lt)(void *data, int i, int j);

/* Must never fail */
typedef void (*heapq_swap)(void *data, int i, int j);

int heapq_heapify(int len, void *heap,
				  heapq_lt lt, heapq_swap swap);

/* Call after replacing the first element in a heapified list */
int heapq_siftup(int len, void *heap, int pos,
				 heapq_lt lt, heapq_swap swap);

/* Must heapify first */
int heapq_heapsort(int len, void *heap,
				   heapq_lt lt, heapq_swap swap);
