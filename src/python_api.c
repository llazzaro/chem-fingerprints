#include <Python.h>

#include "chemfp.h"

static PyObject *
version(PyObject *self, PyObject *args) {
  return PyString_FromString(chemfp_version());
}


/* Slightly renamed so it won't share the same name as strerror(3) */
static PyObject *
strerror_(PyObject *self, PyObject *args) {
  int err;
  if (!PyArg_ParseTuple(args, "i:strerror", &err))
    return NULL;
  return PyString_FromString(chemfp_strerror(err));
}

/*************** Hex fingerprint operations  *************/

static PyObject *
hex_isvalid(PyObject *self, PyObject *args) {
  char *s;
  int len;
  if (!PyArg_ParseTuple(args, "s#:hex_isvalid", &s, &len))
    return NULL;
  return PyInt_FromLong(chemfp_hex_isvalid(len, s));
}

static PyObject *
hex_popcount(PyObject *self, PyObject *args) {
  char *s;
  int len;
  if (!PyArg_ParseTuple(args, "s#:hex_popcount", &s, &len))
    return NULL;
  return PyInt_FromLong(chemfp_hex_popcount(len, s));
}

static PyObject *
hex_intersect_popcount(PyObject *self, PyObject *args) {
  char *s1, *s2;
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
  char *s1, *s2;
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
  char *s1, *s2;
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

/********* Byte fingerprint operations  *************/

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

/*************** Internal validation routines  *************/

static int
bad_num_bits(int num_bits) {
  if (num_bits <= 0) {
    PyErr_SetString(PyExc_TypeError, "num_bits must be positive");
    return 1;
  }
  return 0;
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
  if (threshold < 0.0 || threshold > 1.0) {
    PyErr_SetString(PyExc_TypeError, "threshold must between 0.0 and 1.0, inclusive");
    return 1;
  }
  return 0;
}

/* The arena num bits and storage size must be compatible */
static int
bad_arena_size(const char *which, int num_bits, int storage_size) {
  char msg[150];
  int fp_size = (num_bits+7) / 8;
  if (storage_size < 0) {
    sprintf(msg, "%sstorage_size must be positive", which);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  if (fp_size > storage_size) {
    sprintf(msg, "num_bits of %d (%d bytes) does not fit into %sstorage_size of %d",
            num_bits, fp_size, which, storage_size);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  return 0;
}

/* There must be enough cells for at least num queries (in an FPS threshold search) */
static int 
bad_fps_cells(int *num_cells, int cells_size, int num_queries) {
  char msg[100];
  *num_cells = cells_size / sizeof(chemfp_tanimoto_cell);
  if (*num_cells < num_queries) {
    sprintf(msg, "%d queries requires at least %d cells, not %d",
            num_queries, num_queries, *num_cells);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  return 0;
}

static int
bad_knearest_search_size(int knearest_search_size) {
  if (knearest_search_size < sizeof(chemfp_fps_knearest_search)) {
    PyErr_SetString(PyExc_TypeError,
		    "Not enough space allocated for a chemfp_fps_knearest_search");
    return 1;
  }
  return 0;
}

/* Check/adjust the start and end positions into an FPS block */
static int
bad_block_limits(int block_size, int *start, int *end) {
  if (*start < 0) {
    PyErr_SetString(PyExc_TypeError, "block start must not be negative");
    return 1;
  }
  if (*end == -1 || *end > block_size) {
    *end = block_size;
  } else if (*end < 0) {
    PyErr_SetString(PyExc_TypeError, "block end must either be -1 or non-negative");
    return 1;
  }

  if (*start > block_size) {
    *start = block_size;
  }
  return 0;
}

/* Check/adjust the start and end positions into an arena */
static int
bad_arena_limits(const char *which, int arena_size, int storage_size, int *start, int *end) {
  char msg[150];
  int max_index;
  if (arena_size % storage_size != 0) {
    sprintf(msg, "%sarena size (%d) is not a multiple of its storage size (%d)",
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
  if (*start > max_index) {  /* I'll later ignore if start is too large */
    *start = max_index;
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

static int
bad_fingerprint_sizes(int num_bits, int query_storage_size, int target_storage_size) {
  return (bad_arena_size("query_", num_bits, query_storage_size) ||
          bad_arena_size("target_", num_bits, target_storage_size));
}

static int
bad_popcount_indicies(const char *which, int check_indicies, int num_bits,
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
            which, popcount_indicies_size);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }

  num_popcounts = popcount_indicies_size / sizeof(int);

  if (num_bits > num_popcounts - 1) {
    sprintf(msg, "%d bits requires at least %d %spopcount indicies, not %d",
            num_bits, num_bits+1, which, num_popcounts);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }

  if (check_indicies) {
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
  }
  return 0;
}



static int
bad_offsets(int num_queries, int size, int start) {
  int num_offsets = size/sizeof(int);

  // There must be enough space for num_queries, plus 1 for the end
  char msg[100];
  if (start < 0) {
    PyErr_SetString(PyExc_TypeError, "result offsets start must not be negative");
    return 1;
  }
  if (start > num_offsets) {
    sprintf(msg, "result offsets start (%d) is too large", start);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  if ((num_queries+1) < (num_offsets - start)) {
    sprintf(msg, "Insuffient space to store %d result offsets", num_queries);
    PyErr_SetString(PyExc_TypeError, msg);
    return 1;
  }
  return 0;
}

static int
bad_cells(int min_row_size,int indicies_size, int scores_size, int *num_cells) {
  int num_indicies = indicies_size / sizeof(int);
  int num_scores = scores_size / sizeof(double);

  if (num_indicies < min_row_size) {
    PyErr_SetString(PyExc_TypeError, "Insufficient space to store indicies for a row");
    return 1;
  }
  if (num_scores < min_row_size) {
    PyErr_SetString(PyExc_TypeError, "Insufficient space to store scores for a row");
    return 1;
  }
  if (num_indicies < num_scores) {
    *num_cells = num_indicies;
  } else {
    *num_cells = num_scores;
  }
  return 0;
}

static int
bad_counts(int count_size, int num_queries) {
  if (count_size / sizeof(int) < num_queries) {
    PyErr_SetString(PyExc_TypeError, "Insufficient space to store all of the counts");
    return 1;
  }
  return 0;
}


/*************** FPS functions  *************/
static int
bad_hex_size(int hex_size) {
  if (hex_size == -1) {
    return 0;
  }
  if (hex_size < 1) {
    PyErr_SetString(PyExc_TypeError, "hex_size must be positive or -1");
    return 1;
  }
  if (hex_size % 2 != 0) {
    PyErr_SetString(PyExc_TypeError, "hex_size must be a multiple of 2");
    return 1;
  }
  return 0;
}

// Is this something I really need? Peering into a block might be better
static PyObject *
fps_line_validate(PyObject *self, PyObject *args) {
  int hex_size, line_size;
  char *line;
  
  if (!PyArg_ParseTuple(args, "is#:fps_line_validate", &hex_size, &line, &line_size))
    return NULL;
  if (bad_hex_size(hex_size))
    return NULL;
  return PyInt_FromLong(chemfp_fps_line_validate(hex_size, line_size, line));
}

/* Extract the binary fingerprint and identifier from the line */

/* This assume only the characters 0-9, A-F and a-f will be used */
static const int _hex_digit_to_value[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  0-15 */
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 16-31 */
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 32-47 */
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, /* 48-63 */
  0,10,11,12,13,14,15, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 64-79 */
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 80-95 */
  0,10,11,12,13,14,15};                           /* 96-102 */

static PyObject *
fps_parse_id_fp(PyObject *self, PyObject *args) {
  int hex_size, line_size, err, i;
  char *line;
  const char *id_start, *id_end;
  PyObject *fp, *retval;
  char *s;
  
  if (!PyArg_ParseTuple(args, "is#:fps_line_validate", &hex_size, &line, &line_size))
    return NULL;
  if (bad_hex_size(hex_size))
    return NULL;

  if (line_size == 0 || line[line_size-1] != '\n') {
    return Py_BuildValue("i(ss)", CHEMFP_MISSING_NEWLINE, NULL, NULL);
  }
  err = chemfp_fps_find_id(hex_size, line, &id_start, &id_end);
  if (err != CHEMFP_OK) {
    return Py_BuildValue("i(ss)", err, NULL, NULL);
  }
  if (hex_size == -1) {
    hex_size = (id_start-line)-1;
  }
  fp = PyString_FromStringAndSize(NULL, hex_size/2);
  if (!fp)
    return NULL;
  s = PyString_AS_STRING(fp);
  for (i=0; i<hex_size; i+=2) {
    *s++ = ((_hex_digit_to_value[(int)line[i]]<<4)+_hex_digit_to_value[(int)line[i+1]]);
  }

  retval = Py_BuildValue("i(s#O)", err, id_start, id_end-id_start, fp);

  /* The "O" added a reference which I need to take away */
  Py_DECREF(fp);

  return retval;
}



/* In Python this is
 (err, num_lines_processed) = fps_tanimoto_count(
     num_bits, query_storage_size, query_arena,
     target_block, target_start, target_end,
     threshold, counts)
*/
static PyObject *
fps_tanimoto_count(PyObject *self, PyObject *args) {
  int num_bits, query_storage_size, query_arena_size, query_start, query_end;
  const unsigned char *query_arena;
  const char *target_block;
  int target_block_size, target_start, target_end;
  double threshold;
  int *counts, counts_size;
  int num_lines_processed = 0;
  int err;

  if (!PyArg_ParseTuple(args, "iit#iit#iidw#",
                        &num_bits,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
                        &target_block, &target_block_size,
                        &target_start, &target_end,
                        &threshold,
			&counts, &counts_size))
    return NULL;

  if (bad_num_bits(num_bits) ||
      bad_arena_size("query_", num_bits, query_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_block_limits(target_block_size, &target_start, &target_end) ||
      bad_threshold(threshold) ||
      bad_counts(counts_size, query_arena_size / query_storage_size)) {
    return NULL;
  }
  if (target_start >= target_end) {
    // start of next byte to process, num lines processed, num cells
    return Py_BuildValue("iiii", CHEMFP_OK, 0);
  }
  err = chemfp_fps_tanimoto_count(
        num_bits, 
        query_storage_size, query_arena, query_start, query_end,
        target_block+target_start, target_end-target_start,
        threshold, counts, &num_lines_processed);

  return Py_BuildValue("ii", err, num_lines_processed);
                       
}

/* In Python this is
 (err, next_start, num_lines_processed, num_cells_processed) = 
     fps_threshold_tanimoto_search(num_bits, query_storage_size, query_arena,
                                   target_block, target_start, target_end,
				   threshold, cells)
*/
static PyObject *
fps_threshold_tanimoto_search(PyObject *self, PyObject *args) {
  int num_bits, query_storage_size, query_arena_size, query_start, query_end;
  const unsigned char *query_arena;
  const char *target_block, *stopped_at;
  int target_block_size, target_start, target_end;
  chemfp_tanimoto_cell *cells;
  double threshold;
  int cells_size;
  int num_lines_processed = 0, num_cells_processed = 0;
  int num_cells, err;

  if (!PyArg_ParseTuple(args, "iit#iit#iidw#",
                        &num_bits,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
                        &target_block, &target_block_size,
                        &target_start, &target_end,
                        &threshold,
                        &cells, &cells_size))
    return NULL;

  if (bad_num_bits(num_bits) ||
      bad_arena_size("query_", num_bits, query_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_block_limits(target_block_size, &target_start, &target_end) ||
      bad_threshold(threshold) ||
      bad_fps_cells(&num_cells, cells_size, query_arena_size / query_storage_size)) {
    return NULL;
  }
  if (target_start >= target_end) {
    // start of next byte to process, num lines processed, num cells
    return Py_BuildValue("iiii", CHEMFP_OK, target_end, 0, 0);
  }
  err = chemfp_fps_threshold_tanimoto_search(
        num_bits, 
        query_storage_size, query_arena, query_start, query_end,
        target_block+target_start, target_end-target_start,
        threshold,
        num_cells, cells,
        &stopped_at, &num_lines_processed, &num_cells_processed);

  return Py_BuildValue("iiii", err, stopped_at - target_block,
                       num_lines_processed, num_cells_processed);
}

static PyObject *
fps_knearest_search_init(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;
  int knearest_search_size, num_bits, query_storage_size;
  unsigned const char *query_arena;
  int query_arena_size, query_start, query_end, k;
  double threshold;
  int err;

  if (!PyArg_ParseTuple(args, "w#iit#iiid",
			&knearest_search, &knearest_search_size,
			&num_bits, &query_storage_size,
			&query_arena, &query_arena_size, &query_start, &query_end,
			&k, &threshold))
    return NULL;
    
  if (bad_knearest_search_size(knearest_search_size) ||
      bad_num_bits(num_bits) ||
      bad_arena_size("query_", num_bits, query_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_k(k) ||
      bad_threshold(threshold)) {
    return NULL;
  }
  err = chemfp_fps_knearest_search_init(
	  knearest_search, num_bits, query_storage_size, 
	  query_arena, query_start, query_end,
	  k, threshold);
  if (err) {
    PyErr_SetString(PyExc_TypeError, chemfp_strerror(err));
    return NULL;
  }
  return Py_BuildValue("");
}

static PyObject *
fps_knearest_search_feed(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;  
  int knearest_search_size;
  const char *target_block;
  int target_block_size, target_start, target_end;
  int err;

  if (!PyArg_ParseTuple(args, "w#t#",
			&knearest_search, &knearest_search_size,
			&target_block, &target_block_size, &target_start, &target_end))
    return NULL;

  if (bad_knearest_search_size(knearest_search_size) ||
      bad_block_limits(target_block_size, &target_start, &target_end))
    return NULL;

  err = chemfp_fps_knearest_search_feed(knearest_search, target_block_size, target_block);
  return PyInt_FromLong(err);
}

static PyObject *
fps_knearest_search_finish(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;  
  int knearest_search_size;
  
  if (!PyArg_ParseTuple(args, "w#",
			&knearest_search, &knearest_search_size))
    return NULL;
  if (bad_knearest_search_size(knearest_search_size))
    return NULL;

  chemfp_fps_knearest_search_finish(knearest_search);
  return Py_BuildValue("");
}


static PyObject *
fps_knearest_search_free(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;  
  int knearest_search_size;
  
  if (!PyArg_ParseTuple(args, "w#",
			&knearest_search, &knearest_search_size))
    return NULL;
  if (bad_knearest_search_size(knearest_search_size))
    return NULL;

  chemfp_fps_knearest_search_free(knearest_search);
  return Py_BuildValue("");
}



/**************** The library-based searches **********/

/* reorder_by_popcount */
static PyObject *
reorder_by_popcount(PyObject *self, PyObject *args) {
  int num_bits;
  int storage_size, start, end;
  unsigned char *arena;
  int arena_size;
  PyObject *py_arena = NULL;
  ChemFPOrderedPopcount *ordering;
  int ordering_size;
  int *popcount_indicies, popcount_indicies_size;

  if (!PyArg_ParseTuple(args, "iit#iiw#w#",
                        &num_bits,
                        &storage_size, &arena, &arena_size,
                        &start, &end,
                        &ordering, &ordering_size,
                        &popcount_indicies, &popcount_indicies_size
                        )) {
    return NULL;
  }

  if (bad_num_bits(num_bits) ||
      bad_arena_limits("", arena_size, storage_size, &start, &end) ||
      bad_popcount_indicies("", 0, num_bits, popcount_indicies_size, NULL)) {
    return NULL;
  }
  if ((ordering_size / sizeof(ChemFPOrderedPopcount)) < (end-start)) {
    PyErr_SetString(PyExc_TypeError, "allocated ordering space is too small");
    return NULL;
  }
  // TODO: compute the counts first. If everything is in order then
  // there's no need to allocate a new arena string.
  if (end <= start) {
    py_arena = PyString_FromStringAndSize("", 0);
  } else {
    py_arena = PyString_FromStringAndSize(NULL, (end-start)*storage_size);
  }
  if (!py_arena)
    goto error;

  chemfp_reorder_by_popcount(num_bits, storage_size,
                             arena, start, end,
                             (unsigned char *) PyString_AS_STRING(py_arena),
                             ordering, popcount_indicies);

  return py_arena;

 error:
  if (py_arena) {
    PyObject_Del(py_arena);
  }
  return NULL;
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
      bad_num_bits(num_bits) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indicies("target ", 1, num_bits, 
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
threshold_tanimoto_arena(PyObject *self, PyObject *args) {
  double threshold;
  int num_bits;
  int query_storage_size, query_arena_size, query_start, query_end;
  unsigned char *query_arena;
  int target_storage_size, target_arena_size, target_start, target_end;
  unsigned char *target_arena;

  int *target_popcount_indicies, target_popcount_indicies_size;
  int num_cells;
  int *result_offsets, result_offsets_size, result_offsets_start;
  int *result_indicies, result_indicies_size, result_scores_size;
  double *result_scores;

  int result;

    
  if (!PyArg_ParseTuple(args, "diit#iiit#iit#w#iw#w#",
                        &threshold,
                        &num_bits,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
                        &target_storage_size, &target_arena, &target_arena_size,
                        &target_start, &target_end,
                        &target_popcount_indicies, &target_popcount_indicies_size,
                        &result_offsets, &result_offsets_size, &result_offsets_start,
                        &result_indicies, &result_indicies_size,
                        &result_scores, &result_scores_size
                        ))
    return NULL;

  if (bad_threshold(threshold) ||
      bad_num_bits(num_bits) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indicies("target ", 1, num_bits,
                            target_popcount_indicies_size, &target_popcount_indicies)) {
    return NULL;
  }

  if (bad_offsets(query_end-query_start, result_offsets_size, result_offsets_start) ||
      bad_cells(target_end-target_start, result_indicies_size,
                result_scores_size, &num_cells)) {
    return NULL;
  }

  result = chemfp_threshold_tanimoto_arena(
        threshold,
        num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        target_popcount_indicies,
        result_offsets+result_offsets_start,
        num_cells, result_indicies, result_scores);

  return PyInt_FromLong(result);
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
  int num_cells;
  int *result_offsets, result_offsets_size, result_offsets_start;
  int *result_indicies, result_indicies_size, result_scores_size;
  double *result_scores;

  int result;

    
  if (!PyArg_ParseTuple(args, "idiit#iiit#iit#w#iw#w#",
                        &k, &threshold,
                        &num_bits,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
                        &target_storage_size, &target_arena, &target_arena_size,
                        &target_start, &target_end,
                        &target_popcount_indicies, &target_popcount_indicies_size,
                        &result_offsets, &result_offsets_size, &result_offsets_start,
                        &result_indicies, &result_indicies_size,
                        &result_scores, &result_scores_size
                        ))
    return NULL;

  if (bad_k(k) ||
      bad_threshold(threshold) ||
      bad_num_bits(num_bits) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indicies("target ", 1, num_bits,
                            target_popcount_indicies_size, &target_popcount_indicies)) {
    return NULL;
  }

  if (bad_offsets(query_end-query_start, result_offsets_size, result_offsets_start) ||
      bad_cells(k, result_indicies_size, result_scores_size, &num_cells)) {
    return NULL;
  }

  result = chemfp_klargest_tanimoto_arena(
        k, threshold,
        num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        target_popcount_indicies,
        result_offsets+result_offsets_start,
        num_cells, result_indicies, result_scores);

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
   "fps_line_validate (TODO: document)"},
  {"fps_parse_id_fp", fps_parse_id_fp, METH_VARARGS,
   "fps_parse_id_fp (TODO: document)"},

  {"fps_threshold_tanimoto_search", fps_threshold_tanimoto_search, METH_VARARGS,
   "fps_threshold_tanimoto_search (TODO: document)"},

  {"fps_tanimoto_count", fps_tanimoto_count, METH_VARARGS,
   "fps_tanimoto_count (TODO: document)"},


  {"fps_knearest_search_init", fps_knearest_search_init, METH_VARARGS,
   "fps_knearest_search_init (TODO: document)"},
  {"fps_knearest_search_feed", fps_knearest_search_feed, METH_VARARGS,
   "fps_knearest_search_feed (TODO: document)"},
  {"fps_knearest_search_finish", fps_knearest_search_finish, METH_VARARGS,
   "fps_knearest_search_finish (TODO: document)"},
  {"fps_knearest_search_free", fps_knearest_search_free, METH_VARARGS,
   "fps_knearest_search_free (TODO: document)"},


#if 0
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
#endif

  {"count_tanimoto_arena", count_tanimoto_arena, METH_VARARGS,
   "count_tanimoto_arena (TODO: document)"},

  {"threshold_tanimoto_arena", threshold_tanimoto_arena, METH_VARARGS,
   "threshold_tanimoto_arena (TODO: document)"},

  {"klargest_tanimoto_arena", klargest_tanimoto_arena, METH_VARARGS,
   "klargest_tanimoto_arena (TODO: document)"},

  {"reorder_by_popcount", reorder_by_popcount, METH_VARARGS,
   "reorder_by_popcount (TODO: document)"},

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
