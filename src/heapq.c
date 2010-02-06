#include "heapq.h"

#ifdef TEST_MAIN
#include <stdio.h>
#endif

/* Implement nlargest from the Python */

int heapq_siftdown(int len, void *heap, int startpos, int pos,
				   heapq_lt lt, heapq_swap swap) {
  int parentpos, cmp;
  /* Follow the path to the root, moving parents down until finding
	 a place newitem fits. */
  while (pos > startpos) {
	parentpos = (pos-1) >> 1;
	cmp = lt(heap, pos, parentpos);
	if (cmp == -1)
	  return -1;
	if (cmp == 0)
	  break;
	swap(heap, pos, parentpos);
	pos = parentpos;
  }
  return 0;
}

int heapq_siftup(int len, void *heap, int pos,
				 heapq_lt lt, heapq_swap swap) {
  int endpos = len;
  int startpos = pos;
  int cmp, rightpos, childpos;

  /* Bubble up the smaller child until hitting a leaf. */
  childpos = 2*pos + 1;
  while (childpos < endpos) {
	/* Set childpos to index of smaller child. */
	rightpos = childpos + 1;
	if (rightpos < endpos) {
	  cmp = lt(heap, childpos, rightpos);
	  if (cmp == -1)
		return -1;
	  if (cmp == 0)
		childpos = rightpos;
	}
	/* Move the smaller child up. */
	swap(heap, pos, childpos);
	pos = childpos;
	childpos = 2*pos + 1;
  }
  /* The item at pos contains the original 'pos' item.  Bubble it back up to
	 its final resting place (by sifting its parents down). */
  return heapq_siftdown(len, heap, startpos, pos, lt, swap);
}

int heapq_heapify(int len, void *heap,
				  heapq_lt lt, heapq_swap swap) {
  int i;
  /* Transform bottom-up.  The largest index there's any point to
	 looking at is the largest with a child index in-range, so must
	 have 2*i + 1 < n, or i < (n-1)/2.  If n is even = 2*j, this is
	 (2*j-1)/2 = j-1/2 so j-1 is the largest, which is n//2 - 1.  If
	 n is odd = 2*j+1, this is (2*j+1-1)/2 = j so j-1 is the largest,
	 and that's again n//2-1.
  */
  for (i=len/2-1; i>=0; i--) {
	if (heapq_siftup(len, heap, i, lt, swap) == -1)
	  return -1;
  }
  return 0;
}

/* http://en.wikipedia.org/wiki/Heapsort */

int heapq_heapsort(int len, void *heap,
				   heapq_lt lt, heapq_swap swap) {
  int end;
  if (len == 0)
	return 0;
  for (end = len-1; end>0; end--) {
	swap(heap, 0, end);
	if (heapq_siftup(end, heap, 0, lt, swap) == -1)
	  return -1;
  }
  return 0;
}

#ifdef TEST_MAIN

int int_lt(void *heap, int left, int right) {
  int *data = (int *) heap;
  return data[left] < data[right];
}

void int_swap(void *heap, int left, int right) {
  int *data = (int *) heap;
  int tmp = data[left];
  data[left] = data[right];
  data[right] = tmp;
}

void add_int(int len, int *data, int new_value) {
  if (new_value < data[0]) 
	return;
  data[0] = new_value;
  heapq_siftup(len, data, 0, int_lt, int_swap);
}
	

void dump(int len, int *data) {
  int i;
  for (i=0; i<len; i++) {
	printf("%d ", data[i]);
  }
  printf("\n");
}

int main() {
  int data[] = {3, 9, 2, 4, 5};
  int new_data[] = {8, 6, 7, 5, 3, 0, 9};
  int len = sizeof(data)/sizeof(int);
  int new_len = sizeof(new_data)/sizeof(int);
  int i;

  heapq_heapify(len, data, int_lt, int_swap);
  dump(len, data);
  for (i=0; i<new_len; i++) {
	add_int(len, data, new_data[i]);
	dump(len, data);
  }
  heapq_heapsort(len, data, int_lt, int_swap);
  dump(len, data);

  return 0;
}
#endif
