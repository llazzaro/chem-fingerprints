/* Low-level heap commands */

/* These are private internal functions used by the rest of the chemfp code */
/* They are not part of the public API */

#include "heapq.h"

/* This code is derived from Python's _heapqmodule.c

Heritage from _heapqmodule.c

   Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010
   Python Software Foundation; All Rights Reserved

   C implementation derived directly from heapq.py in Py2.3 which was
   written by Kevin O'Connor, augmented by Tim Peters, annotated by
   François Pinard, and converted to C by Raymond Hettinger.

I'm only using a few of the functions from that module. For those, I stripped
out the dependencies on Python's data structures and changed it to take
user-defined comparison and swap functions. Those are gloss; the code hasn't
really changed at all.

*/

int chemfp_heapq_siftdown(int len, void *heap, int startpos, int pos,
                          chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
  int parentpos, cmp;
  /* unused parameter */
  (void)(len);
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

/* Put the heap into sorted order. The code must already be heapified. */
/* Details at http://en.wikipedia.org/wiki/Heapsort */

int chemfp_heapq_heapsort(int len, void *heap,
                          chemfp_heapq_lt lt, chemfp_heapq_swap swap) {
  int end;
  if (len == 0)
    return 0;
  for (end = len-1; end>0; end--) {
    swap(heap, 0, end);
    if (chemfp_heapq_siftup(end, heap, 0, lt, swap) == -1)
      return -1;
  }
  return 0;
}
