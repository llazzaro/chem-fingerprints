#include <Python.h>

#include "chemfp.h"

static PyObject *
version(PyObject *self, PyObject *args) {
  return PyString_FromString(chemfp_version());
}


// Slightly renamed so it won't share the same name as strerror(3)
static PyObject *
strerror_(PyObject *self, PyObject *args) {
  int err;
  if (!PyArg_ParseTuple(args, "i:strerror", &err))
    return NULL;
  return PyString_FromString(chemfp_strerror(err));
}


static PyObject *
hex_isvalid(PyObject *self, PyObject *args) {
  unsigned char *s;
  int len;
  if (!PyArg_ParseTuple(args, "s#:hex_isvalid", &s, &len))
    return NULL;
  return PyInt_FromLong(chemfp_hex_isvalid(len, s));
}

static PyObject *
hex_popcount(PyObject *self, PyObject *args) {
  unsigned char *s;
  int len;
  if (!PyArg_ParseTuple(args, "s#:hex_popcount", &s, &len))
    return NULL;
  return PyInt_FromLong(chemfp_hex_popcount(len, s));
}

static PyObject *
hex_intersect_popcount(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:hex_intersect_popcount", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "hex fingerprints must have the same length");
    return NULL;
  }
  return PyInt_FromLong(chemfp_hex_intersect_popcount(len1, s1, s2));
}

static PyObject *
hex_tanimoto(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:hex_tanimoto", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "hex fingerprints must have the same length");
    return NULL;
  }
  return PyFloat_FromDouble(chemfp_hex_tanimoto(len1, s1, s2));
}

static PyObject *
hex_contains(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:hex_contains", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "hex fingerprints must have the same length");
    return NULL;
  }
  return PyInt_FromLong(chemfp_hex_contains(len1, s1, s2));
}
/*** ***/

static PyObject *
byte_popcount(PyObject *self, PyObject *args) {
  unsigned char *s;
  int len;
  if (!PyArg_ParseTuple(args, "s#:byte_popcount", &s, &len))
    return NULL;
  return PyInt_FromLong(chemfp_byte_popcount(len, s));
}

static PyObject *
byte_intersect_popcount(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:byte_intersect_popcount", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "byte fingerprints must have the same length");
    return NULL;
  }
  return PyInt_FromLong(chemfp_byte_intersect_popcount(len1, s1, s2));
}

static PyObject *
byte_tanimoto(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:byte_tanimoto", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "byte fingerprints must have the same length");
    return NULL;
  }
  return PyFloat_FromDouble(chemfp_byte_tanimoto(len1, s1, s2));
}

static PyObject *
byte_contains(PyObject *self, PyObject *args) {
  unsigned char *s1, *s2;
  int len1, len2;
  if (!PyArg_ParseTuple(args, "s#s#:byte_contains", &s1, &len1, &s2, &len2))
    return NULL;
  if (len1 != len2) {
    PyErr_SetString(PyExc_TypeError,
                    "byte fingerprints must have the same length");
    return NULL;
  }
  return PyInt_FromLong(chemfp_byte_contains(len1, s1, s2));
}

//

static PyObject *
fps_line_validate(PyObject *self, PyObject *args) {
  int hex_len, line_len;
  char *line;
  
  if (!PyArg_ParseTuple(args, "is#:fps_line_validate", &hex_len, &line, &line_len))
    return NULL;
  return PyInt_FromLong(chemfp_fps_line_validate(hex_len, line_len, line));
}

static PyObject *
fps_tanimoto(PyObject *self, PyObject *args) {
  int hex_len, target_block_len, *id_lens;
  char *hex_query, *target_block, **id_starts;
  double threshold, *scores;
  int *lineno_p = NULL, *num_found;
  int err;
  
  if (!PyArg_ParseTuple(args, "s#s#dwwww|w:fps_tanimoto",
                        &hex_query, &hex_len,
                        &target_block, &target_block_len,
                        &threshold,
                        &num_found,
                        &id_starts, &id_lens, &scores,
                        &lineno_p))
    return NULL;
  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_tanimoto(hex_len, hex_query,
                            target_block_len, target_block,
                            threshold,
                            num_found,
                            id_starts, id_lens, scores,
                            lineno_p);
  Py_END_ALLOW_THREADS;

  return PyInt_FromLong(err);
}

static PyObject *
fps_tanimoto_count(PyObject *self, PyObject *args) {
  int hex_len, target_block_len, *num_found, *lineno=NULL;
  char *hex_query, *target_block;
  double threshold;
  int err;
  if (!PyArg_ParseTuple(args, "s#s#dw|w:fps_tanimoto_count",
                        &hex_query, &hex_len,
                        &target_block, &target_block_len,
                        &threshold,
                        &num_found, &lineno))
    return NULL;
  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_tanimoto_count(hex_len, hex_query,
                                  target_block_len, target_block,
                                  threshold,
                                  num_found, lineno);
  Py_END_ALLOW_THREADS;
  return PyInt_FromLong(err);
}

static PyObject *
fps_heap_init(PyObject *self, PyObject *args) {
  chemfp_heap *heap;
  int k;
  double threshold;
  double *scores;
  int *indicies, *id_lens;
  char **id_starts;
  if (!PyArg_ParseTuple(args, "widwwww:fps_heap_init", &heap,
                        &k, &threshold, &indicies, &scores,
                        &id_starts, &id_lens))
    return NULL;
  chemfp_fps_heap_init(heap, k, threshold, indicies, scores,
                         id_starts, id_lens);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
fps_heap_update_tanimoto(PyObject *self, PyObject *args) {
  chemfp_heap *heap;
  int hex_len, target_block_len;
  char *hex_query, *target_block;
  int err;
  if (!PyArg_ParseTuple(args, "ws#s#",
                        &heap, &hex_query, &hex_len,
                        &target_block, &target_block_len))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_heap_update_tanimoto(heap, hex_len, hex_query,
                                        target_block_len, target_block,
                                        NULL);
  Py_END_ALLOW_THREADS;
  return PyInt_FromLong(err);
}

static PyObject *
fps_heap_finish_tanimoto(PyObject *self, PyObject *args) {
  chemfp_heap *heap;
  if (!PyArg_ParseTuple(args, "w", &heap))
    return NULL;
  chemfp_fps_heap_finish_tanimoto(heap);
  Py_INCREF(Py_None);
  return Py_None;
}


////////////////////////////////////////

static PyObject *
nlargest_tanimoto_block(PyObject *self, PyObject *args) {
  unsigned char *query_fp, *target_block, *indicies_buffer, *scores_buffer;
  int query_len, target_block_len, num_targets, offset, storage_len;
  int indicies_len, scores_len;
  int n;
  double threshold;
  if (!PyArg_ParseTuple(args, "is#s#iidw#w#",
                        &n,
                        &query_fp, &query_len,
                        &target_block, &target_block_len,
                        &offset, &storage_len,
                        &threshold,
                        &indicies_buffer, &indicies_len,
                        &scores_buffer, &scores_len))
    return NULL;

  if (offset < 0) {
    PyErr_SetString(PyExc_TypeError, "offset cannot be negative");
    return NULL;
  }
  if (storage_len < 1) {
    PyErr_SetString(PyExc_TypeError, "storage_len must be positive");
    return NULL;
  }
  if (! (0<= threshold && threshold <= 1.0)) {
    PyErr_SetString(PyExc_TypeError, "threshold must be between 0.0 and 1.0 inclusive");
    return NULL;
  }
  if (n < 1) {
    PyErr_SetString(PyExc_TypeError, "n must be positive");
    return NULL;
  }
  if (query_len > storage_len) {
    PyErr_SetString(PyExc_TypeError,
                    "query fingerprint is longer than target fingerprint storage_len");
    return NULL;
  }

  /* Work out how may fingerprints there are */

  if (offset > target_block_len) {
    PyErr_SetString(PyExc_TypeError,
                    "offset is larger than the target_block buffer");
    return NULL;
  }
  target_block_len -= offset;
  target_block += offset;
  offset = 0;
  if (target_block_len == 0)
    return PyInt_FromLong(0);

  if ((target_block_len % storage_len) != 0) {
    PyErr_SetString(PyExc_TypeError,
                    "adjusted target_block length is not a multiple of the storage size");
    return NULL;
  }
  num_targets = target_block_len / storage_len;
  if (n > num_targets)
    n = num_targets;

  if (n * sizeof(int) > indicies_len) {
    PyErr_SetString(PyExc_TypeError, "indicies buffer is not long enough");
    return NULL;
  }
  if (n * sizeof(double) > scores_len) {
    PyErr_SetString(PyExc_TypeError, "score buffer is not long enough");
    return NULL;
  }

  return PyInt_FromLong(chemfp_nlargest_tanimoto_block(
        n,
        query_len, query_fp,
        num_targets, target_block, offset, storage_len,
        threshold,
        (int *) indicies_buffer, (double *) scores_buffer));
}

static PyObject *
hex_nlargest_tanimoto_block(PyObject *self, PyObject *args) {
  int n, query_len, target_block_len;
  unsigned char *query_fp, *target_block;
  double threshold;
  int endpos;
  double *scores;
  int *id_lens, *lineno;
  unsigned char **start_ids;
  int err;
  if (!PyArg_ParseTuple(args, "is#s#idwwww:hex_largest_tanimoto_block",
                        &n, &query_fp, &query_len,
                        &target_block, &target_block_len,
                        &endpos,
                        &threshold,
                        &scores,
                        &start_ids,
                        &id_lens,
                        &lineno))
    return NULL;

  err = chemfp_hex_tanimoto_block(n, query_len, query_fp,
                                  endpos, target_block,
                                  threshold,
                                  scores, start_ids, id_lens, lineno);
  return PyInt_FromLong(err);

}
static PyObject *
intersect_popcount_count(PyObject *self, PyObject *args) {
  unsigned char *query_fp, *target_block;
  int query_len, target_block_len, offset, storage_len, min_overlap;
  int num_targets;
  if (!PyArg_ParseTuple(args, "s#s#iii",
						&query_fp, &query_len,
						&target_block, &target_block_len,
						&offset, &storage_len,
						&min_overlap))
	return NULL;
  if (offset < 0) {
    PyErr_SetString(PyExc_TypeError, "offset cannot be negative");
    return NULL;
  }
  if (storage_len < 1) {
    PyErr_SetString(PyExc_TypeError, "storage_len must be positive");
    return NULL;
  }
  if (min_overlap < 0) {
	PyErr_SetString(PyExc_TypeError, "min_overlap must be non-negative");
  }
  if (query_len > storage_len) {
    PyErr_SetString(PyExc_TypeError,
                    "query fingerprint is longer than target fingerprint storage_len");
    return NULL;
  }
  if (offset > target_block_len) {
    PyErr_SetString(PyExc_TypeError,
                    "offset is larger than the target_block buffer");
    return NULL;
  }
  target_block_len -= offset;
  target_block += offset;
  offset = 0;
  if (target_block_len == 0)
    return PyInt_FromLong(0);

  if ((target_block_len % storage_len) != 0) {
    PyErr_SetString(PyExc_TypeError,
                    "adjusted target_block length is not a multiple of the storage size");
    return NULL;
  }
  num_targets = target_block_len / storage_len;

  return PyInt_FromLong(
		chemfp_byte_intersect_popcount_count(
		  query_len, query_fp,
          num_targets, target_block, offset, storage_len,
		  min_overlap)
					   );
}

static int
bad_k(int k) {
  if (k < 0) {
    PyErr_SetString(PyExc_TypeError, "k must not be negative");
    return 1;
  }
  return 0;
}

static int
bad_threshold(double threshold) {
  if (threshold <= 0.0 || threshold > 1.0) {
    PyErr_SetString(PyExc_TypeError, "threshold must between 0.0 and 1.0, inclusive");
    return 1;
  }
  return 0;
}

static int
bad_fingerprint_sizes(int num_bits, int query_storage_size, int target_storage_size) {
  char msg[150];
  int fp_size = num_bits / 8;
  if (num_bits <= 0) {
    PyErr_SetString(PyExc_TypeError, "num_bits must be positive");
    return 1;
  }
  if (query_storage_size < 0) {
    PyErr_SetString(PyExc_TypeError, "query_storage_size must be positive");
    return 1;
  }
  if (target_storage_size < 0) {
    PyErr_SetString(PyExc_TypeError, "target_storage_size must be positive");
    return 1;
  }

  if (fp_size > query_storage_size) {
    sprintf(msg, "num_bits of %d (%d bytes) does not fit into query_storage_size of %d",
	    num_bits, fp_size, query_storage_size);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  if (fp_size > target_storage_size) {
    sprintf(msg, "num_bits of %d (%d bytes) does not fit into target_storage_size of %d",
	    num_bits, fp_size, target_storage_size);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  return 0;
}

static int
bad_popcount_indicies(const char *which, int num_bits,
		      int popcount_indicies_size, int **popcount_indicies_ptr) {
  char msg[150];
  int num_popcounts;
  int prev, i;
  int *popcount_indicies;

  if (popcount_indicies_size == 0) {
    /* Special case: this means to ignore this field */
    *popcount_indicies_ptr = NULL;
    return 0;
  }
  if ((popcount_indicies_size % sizeof(int)) != 0) {
    sprintf(msg,
	    "%spopcount indicies length (%d) is not a multiple of the native integer size",
	    which, num_bits);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }

  num_popcounts = popcount_indicies_size / sizeof(int);

  if (num_bits < num_popcounts - 1) {
    sprintf(msg, "%d bits requires at least %d %spopcount indicies, not %d",
	    num_bits, num_bits+1, which, num_popcounts);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }

  popcount_indicies = *popcount_indicies_ptr;
  if (popcount_indicies[0] != 0) {
    sprintf(msg, "%s popcount indicies[0] must be 0", which);
    PyErr_SetString(PyExc_TypeError, "%spopcount_indicies[0] must be 0");
      return 1;
  }
  prev = 0;
  for (i=1; i<num_popcounts; i++) {
    if (popcount_indicies[i] < prev) {
      sprintf(msg, "%spopcount indicies must never decrease", which);
      PyErr_SetString(PyExc_TypeError, msg);
      return 1;
    }
    prev = popcount_indicies[i];
  }
  return 0;
}

static int
bad_limits(const char *which, int arena_size, int storage_size, int *start, int *end) {
  char msg[150];
  int max_index;
  if (arena_size % storage_size != 0) {
    sprintf(msg, "%s arena size (%d) is not a multiple of its storage size (%d)",
	    which, arena_size, storage_size);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  if (*start < 0) {
    sprintf(msg, "%sstart must not be negative", which);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  max_index = arena_size / storage_size;
  if (*start > max_index) {
    *start = max_index; // XXX Why isn't this an error?
  }
  if (*end == -1 || *end > max_index) {
    *end = max_index;
  } else if (*end < 0) {
    sprintf(msg, "%send must either be -1 or non-negative", which);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  return 0;
}



/* count_tanimoto_arena */
static PyObject *
count_tanimoto_arena(PyObject *self, PyObject *args) {
  double threshold;
  int num_bits;
  unsigned char *query_arena, *target_arena;
  int query_storage_size, query_arena_size=0, query_start=0, query_end=0;
  int target_storage_size, target_arena_size=0, target_start=0, target_end=0;
  
  int *target_popcount_indicies, target_popcount_indicies_size;
  
  int result_counts_size, *result_counts;

  if (!PyArg_ParseTuple(args, "diis#iiis#iis#w#",
			&threshold,
			&num_bits,
			&query_storage_size, &query_arena, &query_arena_size,
			&query_start, &query_end,
			&target_storage_size, &target_arena, &target_arena_size,
			&target_start, &target_end,
			&target_popcount_indicies, &target_popcount_indicies_size,
			&result_counts, &result_counts_size))
    return NULL;

  if (bad_threshold(threshold) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_limits("query ", query_arena_size, query_storage_size,
		 &query_start, &query_end) ||
      bad_limits("target ", target_arena_size, target_storage_size,
		 &target_start, &target_end) ||
      bad_popcount_indicies("target ", num_bits,
			    target_popcount_indicies_size, &target_popcount_indicies)) {
    return NULL;
  }

  if (query_start > query_end) {
    Py_RETURN_NONE;
  }

  if (result_counts_size < (query_end - query_start)*sizeof(int)) {
    PyErr_SetString(PyExc_TypeError, "not enough space allocated for result_counts");
    return NULL;
  }

  chemfp_count_tanimoto_arena(threshold,
			      num_bits,
			      query_storage_size, query_arena, query_start, query_end,
			      target_storage_size, target_arena, target_start, target_end,
			      target_popcount_indicies,
			      result_counts);
  Py_RETURN_NONE;
}
    

/* klargest_tanimoto_arena */
static PyObject *
klargest_tanimoto_arena(PyObject *self, PyObject *args) {
  int k;
  double threshold;
  int num_bits;
  int query_storage_size, query_arena_size, query_start, query_end;
  unsigned char *query_arena;
  int target_storage_size, target_arena_size, target_start, target_end;
  unsigned char *target_arena;

  int *target_popcount_indicies, target_popcount_indicies_size;
  int num_allocated;
  int *result_counts, *result_indicies;
  double *result_scores;

  int result;

    
  if (!PyArg_ParseTuple(args, "idiis#iis#iit#iwww",
			&k, &threshold,
			&num_bits,
			&query_storage_size, &query_arena, &query_arena_size,
			&query_start, &query_end,
			&target_storage_size, &target_arena, &target_arena_size,
			&target_start, &target_end,
			&target_popcount_indicies, &target_popcount_indicies_size,
			&num_allocated,
			&result_counts, &result_indicies, &result_scores
			))
    return NULL;

  if (bad_k(k) ||
      bad_threshold(threshold) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_limits("query ", query_arena_size, query_storage_size,
		 &query_start, &query_end) ||
      bad_limits("target ", target_arena_size, target_storage_size,
		 &query_start, &query_end) ||
      bad_popcount_indicies("target ", num_bits,
			    target_popcount_indicies_size, &target_popcount_indicies)) {
    return NULL;
  }
      
  if (num_allocated < 0) {
    PyErr_SetString(PyExc_TypeError, "num_allocated must not be negative");
    return NULL;
  }

  result = chemfp_klargest_tanimoto_arena(
	k, threshold,
	num_bits,
	query_storage_size, query_arena, query_start, query_end,
	target_storage_size, target_arena, target_start, target_end,
	target_popcount_indicies,
	num_allocated, result_counts, result_indicies, result_scores);

  return PyInt_FromLong(result);
}


static PyMethodDef chemfp_methods[] = {
  {"version", version, METH_NOARGS,
   "version()\n\nReturn the chemfp library version, as a string like '1.0'"},
  {"strerror", strerror_, METH_VARARGS,
   "strerror(n)\n\nConvert the error code integer to more descriptive text"},
  {"hex_isvalid", hex_isvalid, METH_VARARGS,
   "hex_isvalid(s)\n\nReturn 1 if the string is a valid hex fingerprint, otherwise 0"},
  {"hex_popcount", hex_popcount, METH_VARARGS, 
   "hex_popcount(fp)\n\nReturn the number of bits set in a hex fingerprint, or -1 for non-hex strings"},
  {"hex_intersect_popcount", hex_intersect_popcount, METH_VARARGS,
   "hex_intersect_popcount(fp1, fp2)\n\nReturn the number of bits set in the intersection of the two hex fingerprint,\nor -1 if either string is a non-hex string"},
  {"hex_tanimoto", hex_tanimoto, METH_VARARGS,
   "hex_tanimoto(fp1, fp2)\n\nCompute the Tanimoto similarity between two hex fingerprints.\nReturn a float between 0.0 and 1.0, or -1.0 if either string is not a hex fingerprint"},
  {"hex_contains", hex_contains, METH_VARARGS,
   "hex_contains(super_fp, sub_fp)\n\nReturn 1 if the on bits of sub_fp are also 1 bits in super_fp, otherwise 0.\nReturn -1 if either string is not a hex fingerprint"},


  {"byte_popcount", byte_popcount, METH_VARARGS,
   "byte_popcount(fp)\n\nReturn the number of bits set in a byte fingerprint"},
  {"byte_intersect_popcount", byte_intersect_popcount, METH_VARARGS,
   "byte_intersect_popcount(fp1, fp2)\n\nReturn the number of bits set in the instersection of the two byte fingerprints"},
  {"byte_tanimoto", byte_tanimoto, METH_VARARGS,
   "byte_tanimoto(fp1, fp2)\n\nCompute the Tanimoto similarity between two byte fingerprints"},
  {"byte_contains", byte_contains, METH_VARARGS,
   "byte_contains(super_fp, sub_fp)\n\nReturn 1 if the on bits of sub_fp are also 1 bits in super_fp"},

  // FPS
  {"fps_line_validate", fps_line_validate, METH_VARARGS,
   "fps_line_validate(s)\n\nReturn 1 if the string is a valid FPS line, else return 0"},
  {"fps_tanimoto", fps_tanimoto, METH_VARARGS,
   "Calculate Tanimoto scores against a block of FPS lines (TODO: document)"},
  {"fps_tanimoto_count", fps_tanimoto_count, METH_VARARGS,
   "Count Tanimoto scores (TODO: document)"},

  {"fps_heap_init", fps_heap_init, METH_VARARGS, "init heap (TODO: document)"},
  {"fps_heap_update_tanimoto", fps_heap_update_tanimoto, METH_VARARGS, 
   "update heap (TODO: document)"},
  {"fps_heap_finish_tanimoto", fps_heap_finish_tanimoto, METH_VARARGS,
   "finish heap (TODO: document)"},

  {"nlargest_tanimoto_block", nlargest_tanimoto_block, METH_VARARGS,
   "nlargest_tanimoto_block (TODO: document)"},

  {"hex_nlargest_tanimoto_block", hex_nlargest_tanimoto_block, METH_VARARGS,
   "hex_nlargest_tanimoto_block (TODO: document)"},

  {"intersect_popcount_count", intersect_popcount_count, METH_VARARGS,
   "intersect_popcount_count (TODO: document)"},

  {"count_tanimoto_arena", count_tanimoto_arena, METH_VARARGS,
   "count_tanimoto_arena (TODO: document)"},

  {"klargest_tanimoto_arena", klargest_tanimoto_arena, METH_VARARGS,
   "klargest_tanimoto_arena (TODO: document)"},

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
