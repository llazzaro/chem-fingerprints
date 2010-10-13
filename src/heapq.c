
/* Implement nlargest from the Python */

int chemfp_heapq_siftdown(int len, void *heap, int startpos, int pos,
						  chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
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

int chemfp_heapq_siftup(int len, void *heap, int pos,
						chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
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
  return chemfp_heapq_siftdown(len, heap, startpos, pos, lt, swap);
}

#if 0

/* Convert an unsorted heapt into a heap. Used during testing. */

int chemfp_heapq_heapify(int len, void *heap,
						 chemfp_heapq_lt lt, chemfp_heapq_swap swap);

int chemfp_heapq_heapify(int len, void *heap,
						 chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
  int i;
  /* Transform bottom-up.  The largest index there's any point to
	 looking at is the largest with a child index in-range, so must
	 have 2*i + 1 < n, or i < (n-1)/2.  If n is even = 2*j, this is
	 (2*j-1)/2 = j-1/2 so j-1 is the largest, which is n//2 - 1.  If
	 n is odd = 2*j+1, this is (2*j+1-1)/2 = j so j-1 is the largest,
	 and that's again n//2-1.
  */
  for (i=len/2-1; i>=0; i--) {
	if (chemfp_heapq_siftup(len, heap, i, lt, swap) == -1)
	  return -1;
  }
  return 0;
}
#endif

/* http://en.wikipedia.org/wiki/Heapsort */

int chemfp_heapq_heapsort(int len, void *heap,
						  chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
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
