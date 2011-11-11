#include <Python.h>

#include "chemfp.h"
#include "chemfp_internal.h"

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
    PyErr_SetString(PyExc_ValueError,
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
    PyErr_SetString(PyExc_ValueError,
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
    PyErr_SetString(PyExc_ValueError,
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
    PyErr_SetString(PyExc_ValueError,
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
    PyErr_SetString(PyExc_ValueError,
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
    PyErr_SetString(PyExc_ValueError,
                    "byte fingerprints must have the same length");
    return NULL;
  }
  return PyInt_FromLong(chemfp_byte_contains(len1, s1, s2));
}

/*************** Internal validation routines  *************/

static int
bad_num_bits(int num_bits) {
  if (num_bits <= 0) {
    PyErr_SetString(PyExc_ValueError, "num_bits must be positive");
    return 1;
  }
  return 0;
}

static int
bad_k(int k) {
  if (k < 0) {
    PyErr_SetString(PyExc_ValueError, "k must not be negative");
    return 1;
  }
  return 0;
}

static int
bad_threshold(double threshold) {
  if (threshold < 0.0 || threshold > 1.0) {
    PyErr_SetString(PyExc_ValueError, "threshold must between 0.0 and 1.0, inclusive");
    return 1;
  }
  return 0;
}

static int
bad_alignment(int alignment) {
  if (chemfp_byte_popcount(sizeof(int), (unsigned char *) &alignment) != 1) {
    PyErr_SetString(PyExc_ValueError, "alignment must be a positive power of two");
    return 1;
  }
  return 0;
}

static int
bad_padding(const char *which, int start_padding, int end_padding,
	    const unsigned char **arena, int *arena_size) {
  char msg[150];
  /*  printf("PADDING: %d %d for %d\n", start_padding, end_padding, *arena_size);*/
  if (start_padding < 0) {
    sprintf(msg, "%sstart_padding must not be negative", which);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  if (end_padding < 0) {
    sprintf(msg, "%send_padding must not be negative", which);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  if ((start_padding + end_padding) > *arena_size) {
    sprintf(msg, "%sarena_size is too small for the paddings", which);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  *arena += start_padding;
  *arena_size -= (start_padding + end_padding);
  return 0;
}



/* The arena num bits and storage size must be compatible */
static int
bad_arena_size(const char *which, int num_bits, int storage_size) {
  char msg[150];
  int fp_size = (num_bits+7) / 8;
  if (storage_size < 0) {
    sprintf(msg, "%sstorage_size must be positive", which);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  if (fp_size > storage_size) {
    sprintf(msg, "num_bits of %d (%d bytes) does not fit into %sstorage_size of %d",
            num_bits, fp_size, which, storage_size);
    PyErr_SetString(PyExc_ValueError, msg);
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
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  return 0;
}

static int
bad_knearest_search_size(int knearest_search_size) {
  if (knearest_search_size < sizeof(chemfp_fps_knearest_search)) {
    PyErr_SetString(PyExc_ValueError,
		    "Not enough space allocated for a chemfp_fps_knearest_search");
    return 1;
  }
  return 0;
}

/* Check/adjust the start and end positions into an FPS block */
static int
bad_block_limits(int block_size, int *start, int *end) {
  if (*start < 0) {
    PyErr_SetString(PyExc_ValueError, "block start must not be negative");
    return 1;
  }
  if (*end == -1 || *end > block_size) {
    *end = block_size;
  } else if (*end < 0) {
    PyErr_SetString(PyExc_ValueError, "block end must either be -1 or non-negative");
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
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  if (*start < 0) {
    sprintf(msg, "%sstart must not be negative", which);
    PyErr_SetString(PyExc_ValueError, msg);
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
    PyErr_SetString(PyExc_ValueError, msg);
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
bad_popcount_indices(const char *which, int check_indices, int num_bits,
                      int popcount_indices_size, int **popcount_indices_ptr) {
  char msg[150];
  int num_popcounts;
  int prev, i;
  int *popcount_indices;

  if (popcount_indices_size == 0) {
    /* Special case: this means to ignore this field */
    *popcount_indices_ptr = NULL;
    return 0;
  }
  if ((popcount_indices_size % sizeof(int)) != 0) {
    sprintf(msg,
            "%spopcount indices length (%d) is not a multiple of the native integer size",
            which, popcount_indices_size);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }

  /* If there is 1 bit then there must be three indices: */
  /*   indices[0]...indices[1] ==> fingerprints with 0 bits set */
  /*   indices[1]...indices[2] ==> fingerprints with 1 bit set */

  num_popcounts = popcount_indices_size / sizeof(int);

  if (num_bits > num_popcounts - 2) {
    sprintf(msg, "%d bits requires at least %d %spopcount indices, not %d",
            num_bits, num_bits+2, which, num_popcounts);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }

  if (check_indices) {
    popcount_indices = *popcount_indices_ptr;
    if (popcount_indices[0] != 0) {
      sprintf(msg, "%s popcount indices[0] must be 0", which);
      PyErr_SetString(PyExc_ValueError, "%spopcount_indices[0] must be 0");
      return 1;
    }
    prev = 0;
    for (i=1; i<num_popcounts; i++) {
      if (popcount_indices[i] < prev) {
        sprintf(msg, "%spopcount indices must never decrease", which);
        PyErr_SetString(PyExc_ValueError, msg);
        return 1;
      }
      prev = popcount_indices[i];
    }
  }
  return 0;
}



static int
bad_offsets(int num_queries, int size, int start) {
  int num_offsets = size/sizeof(int);

  /* There must be enough space for num_queries, plus 1 for the end */
  char msg[100];
  if (start < 0) {
    PyErr_SetString(PyExc_ValueError, "result offsets start must not be negative");
    return 1;
  }
  if (start > num_offsets) {
    sprintf(msg, "result offsets start (%d) is too large", start);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  if ((num_queries+1) < (num_offsets - start)) {
    sprintf(msg, "Insuffient space to store %d result offsets", num_queries);
    PyErr_SetString(PyExc_ValueError, msg);
    return 1;
  }
  return 0;
}

static int
bad_cells(int min_row_size,int indices_size, int scores_size, int *num_cells) {
  int num_indices = indices_size / sizeof(int);
  int num_scores = scores_size / sizeof(double);

  if (num_indices < min_row_size) {
    PyErr_SetString(PyExc_ValueError, "Insufficient space to store indices for a row");
    return 1;
  }
  if (num_scores < min_row_size) {
    PyErr_SetString(PyExc_ValueError, "Insufficient space to store scores for a row");
    return 1;
  }
  if (num_indices < num_scores) {
    *num_cells = num_indices;
  } else {
    *num_cells = num_scores;
  }
  return 0;
}

static int
bad_counts(int count_size, int num_queries) {
  if (count_size / sizeof(int) < num_queries) {
    PyErr_SetString(PyExc_ValueError, "Insufficient space to store all of the counts");
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
    PyErr_SetString(PyExc_ValueError, "hex_size must be positive or -1");
    return 1;
  }
  if (hex_size % 2 != 0) {
    PyErr_SetString(PyExc_ValueError, "hex_size must be a multiple of 2");
    return 1;
  }
  return 0;
}

/* Is this something I really need? Peering into a block might be better */
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
fps_count_tanimoto_hits(PyObject *self, PyObject *args) {
  int num_bits, query_storage_size, query_arena_size, query_start, query_end;
  int query_start_padding, query_end_padding;
  const unsigned char *query_arena;
  const char *target_block;
  int target_block_size, target_start, target_end;
  double threshold;
  int *counts, counts_size;
  int num_lines_processed = 0;
  int err;

  if (!PyArg_ParseTuple(args, "iiiit#iit#iidw#",
                        &num_bits,
			&query_start_padding, &query_end_padding,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
                        &target_block, &target_block_size,
                        &target_start, &target_end,
                        &threshold,
			&counts, &counts_size))
    return NULL;

  if (bad_num_bits(num_bits) ||
      bad_padding("query_", query_start_padding, query_end_padding,
		  &query_arena, &query_arena_size) ||
      bad_arena_size("query_", num_bits, query_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_block_limits(target_block_size, &target_start, &target_end) ||
      bad_threshold(threshold) ||
      bad_counts(counts_size, query_arena_size / query_storage_size)) {
    return NULL;
  }

  if (target_start >= target_end) {
    /* start of next byte to process, num lines processed, num cells */
    return Py_BuildValue("iiii", CHEMFP_OK, 0);
  }
  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_count_tanimoto_hits(
        num_bits, 
        query_storage_size, query_arena, query_start, query_end,
        target_block+target_start, target_end-target_start,
        threshold, counts, &num_lines_processed);
  Py_END_ALLOW_THREADS;

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
    /* start of next byte to process, num lines processed, num cells */
    return Py_BuildValue("iiii", CHEMFP_OK, target_end, 0, 0);
  }
  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_threshold_tanimoto_search(
        num_bits, 
        query_storage_size, query_arena, query_start, query_end,
        target_block+target_start, target_end-target_start,
        threshold,
        num_cells, cells,
        &stopped_at, &num_lines_processed, &num_cells_processed);
  Py_END_ALLOW_THREADS;

  return Py_BuildValue("iiii", err, stopped_at - target_block,
                       num_lines_processed, num_cells_processed);
}

static PyObject *
fps_knearest_search_init(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;
  int start_padding, end_padding;
  int knearest_search_size, num_bits, query_storage_size;
  unsigned const char *query_arena;
  int query_arena_size, query_start, query_end, k;
  double threshold;
  int err;

  if (!PyArg_ParseTuple(args, "w#iiiit#iiid",
			&knearest_search, &knearest_search_size,
			&num_bits, &start_padding, &end_padding, &query_storage_size,
			&query_arena, &query_arena_size, &query_start, &query_end,
			&k, &threshold))
    return NULL;

  if (bad_knearest_search_size(knearest_search_size) ||
      bad_num_bits(num_bits) ||
      bad_padding("", start_padding, end_padding, &query_arena, &query_arena_size) ||
      bad_arena_size("query_", num_bits, query_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_k(k) ||
      bad_threshold(threshold)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_knearest_search_init(
	  knearest_search, num_bits, query_storage_size, 
	  query_arena, query_start, query_end,
	  k, threshold);
  Py_END_ALLOW_THREADS;
  if (err) {
    PyErr_SetString(PyExc_ValueError, chemfp_strerror(err));
    return NULL;
  }
  return Py_BuildValue("");
}

static PyObject *
fps_knearest_tanimoto_search_feed(PyObject *self, PyObject *args) {
  chemfp_fps_knearest_search *knearest_search;  
  int knearest_search_size;
  const char *target_block;
  int target_block_size, target_start, target_end;
  int err;

  if (!PyArg_ParseTuple(args, "w#t#ii",
			&knearest_search, &knearest_search_size,
			&target_block, &target_block_size, &target_start, &target_end))
    return NULL;

  if (bad_knearest_search_size(knearest_search_size) ||
      bad_block_limits(target_block_size, &target_start, &target_end))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  err = chemfp_fps_knearest_tanimoto_search_feed(knearest_search, target_block_size, target_block);
  Py_END_ALLOW_THREADS;
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

  Py_BEGIN_ALLOW_THREADS;
  chemfp_fps_knearest_search_finish(knearest_search);
  Py_END_ALLOW_THREADS;

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

  Py_BEGIN_ALLOW_THREADS;
  chemfp_fps_knearest_search_free(knearest_search);
  Py_END_ALLOW_THREADS;

  return Py_BuildValue("");
}



/**************** The library-based searches **********/

/* Always allocate space */
static PyObject *
_alloc_aligned_arena(Py_ssize_t size, int alignment,
		     int *start_padding, int *end_padding) {
  PyObject *new_py_string;
  char *s;
  int i;

  new_py_string = PyString_FromStringAndSize(NULL, size+alignment-1);
  if (!new_py_string) {
    return NULL;
  }
  s = PyString_AS_STRING(new_py_string);
  i = ALIGNMENT(s, alignment);
  if (i == 0) {
    *start_padding = 0;
    *end_padding = alignment-1;
  } else {
    *start_padding = alignment - i;
    *end_padding = i-1;
  }
  memset(s, 0, *start_padding);
  memset(s+size+*start_padding, 0, *end_padding);
  return new_py_string;
}

static PyObject *
_align_arena(PyObject *input_arena_obj, int alignment,
	     int *start_padding, int *end_padding) {
  const char *input_arena;
  char *output_arena;
  Py_ssize_t input_arena_size;
  int i;
  PyObject *output_arena_obj;

  if (PyObject_AsCharBuffer(input_arena_obj, &input_arena, &input_arena_size)) {
    PyErr_SetString(PyExc_ValueError, "arena must be a character buffer");
    return NULL;
  }
  i = ALIGNMENT(input_arena, alignment);
  
  /* Already aligned */
  if (i == 0) {
    *start_padding = 0;
    *end_padding = 0;
    Py_INCREF(input_arena_obj);
    return input_arena_obj;
  }

  /* Not aligned. We'll have to move it to a new string */
  output_arena_obj = _alloc_aligned_arena(input_arena_size, alignment,
					  start_padding, end_padding);
  output_arena = PyString_AS_STRING(output_arena_obj);

  /* Copy over into the new string */
  memcpy(output_arena+*start_padding, input_arena, input_arena_size);

  return output_arena_obj;
}

static PyObject *
make_unsorted_aligned_arena(PyObject *self, PyObject *args) {
  int alignment;
  int start_padding, end_padding;
  PyObject *input_arena_obj, *output_arena_obj;

  if (!PyArg_ParseTuple(args, "Oi", &input_arena_obj, &alignment)) {
    return NULL;
  }
  if (bad_alignment(alignment)) {
    return NULL;
  }
  output_arena_obj = _align_arena(input_arena_obj, alignment,
				  &start_padding, &end_padding);
  if (!output_arena_obj) {
    return NULL;
  }
  return Py_BuildValue("iiO", start_padding, end_padding, output_arena_obj);
}


static int
calculate_arena_popcounts(int num_bits, int storage_size, const unsigned char *arena,
			  int num_fingerprints, ChemFPOrderedPopcount *ordering) {
  chemfp_popcount_f calc_popcount;
  const unsigned char *fp;
  int fp_index, popcount, prev_popcount;
  /* Compute the popcounts. (Alignment isn't that important here.) */

  calc_popcount = chemfp_select_popcount(num_bits, storage_size, arena);
  fp = arena;
  for (fp_index = 0; fp_index < num_fingerprints; fp_index++, fp += storage_size) {
    popcount = calc_popcount(storage_size, fp);
    ordering[fp_index].popcount = popcount;
    ordering[fp_index].index = fp_index;
  }

  /* Check if the values are already ordered */

  prev_popcount = ordering[0].popcount;
  for (fp_index = 1; fp_index < num_fingerprints; fp_index++) {
    if (ordering[fp_index].popcount < prev_popcount) {
      return 1; /* Need to sort */
    }
    prev_popcount = ordering[fp_index].popcount;
  }
  return 0; /* Don't need to sort */
}


static int compare_by_popcount(const void *left_p, const void *right_p) {
  const ChemFPOrderedPopcount *left = (ChemFPOrderedPopcount *) left_p;
  const ChemFPOrderedPopcount *right = (ChemFPOrderedPopcount *) right_p;
  if (left->popcount < right->popcount) {
    return -1;
  }
  if (left->popcount > right->popcount) {
    return 1;
  }
  if (left->index < right->index) {
    return -1;
  }
  if (left->index > right->index) {
    return 1;
  }
  return 0;
}


static void
set_popcount_indicies(int num_fingerprints, int num_bits,
		      ChemFPOrderedPopcount *ordering, int *popcount_indices) {
  int popcount, i;

  /* We've sorted by popcount so this isn't so difficult */
  popcount = 0;
  popcount_indices[0] = 0;
  for (i=0; i<num_fingerprints; i++) {
    while (popcount < ordering[i].popcount) {
      popcount++;
      popcount_indices[popcount] = i;
      if (popcount == num_bits) {
	/* We are at or above the limit. We can stop now. */
	i = num_fingerprints;
	break;
	/* Note: with corrupted data it is possible
	   that ->popcount can be > num_bits. This is
	   undefined behavior. I get to do what I want.
	   I decided to treat them as having "max_popcount" bits.
	   After all, I don't want corrupt data to crash the
	   system, and no one is going to validate the input
	   fingerprints for correctness each time.  */
      }
    }
  }
  /* Finish up the high end */
  while (popcount <= num_bits) {
    popcount_indices[++popcount] = num_fingerprints;
  }
}


static PyObject *
make_sorted_aligned_arena(PyObject *self, PyObject *args) {
  int start = 0;
  int num_bits, storage_size, num_fingerprints, ordering_size, popcount_indices_size;
  int start_padding, end_padding;
  PyObject *input_arena_obj, *output_arena_obj;
  const unsigned char *input_arena;
  unsigned char *output_arena;
  Py_ssize_t input_arena_size;
  ChemFPOrderedPopcount *ordering;
  int *popcount_indices;
  int need_to_sort, i;
  int alignment;
  
  if (!PyArg_ParseTuple(args, "iiOiw#w#i",
                        &num_bits,
                        &storage_size, &input_arena_obj,
                        &num_fingerprints,
                        &ordering, &ordering_size,
                        &popcount_indices, &popcount_indices_size,
			&alignment
                        )) {
    return NULL;
  }

  if (PyObject_AsCharBuffer(input_arena_obj,
			     (const char **) &input_arena, &input_arena_size)) {
    PyErr_SetString(PyExc_ValueError, "arena must be a character buffer");
    return NULL;
  }
  if (bad_num_bits(num_bits) ||
      bad_arena_limits("", input_arena_size, storage_size, &start, &num_fingerprints) ||
      bad_popcount_indices("", 0, num_bits, popcount_indices_size, NULL)) {
    return NULL;
  }
  if ((ordering_size / sizeof(ChemFPOrderedPopcount)) < num_fingerprints) {
    PyErr_SetString(PyExc_ValueError, "allocated ordering space is too small");
    return NULL;
  }

  /* Handle the trivial case of no fingerprints */

  if (num_fingerprints == 0) {
    return Py_BuildValue("iiO", 0, 0, input_arena_obj);
  }


  need_to_sort = calculate_arena_popcounts(num_bits, storage_size, input_arena,
					   num_fingerprints, ordering);

  if (!need_to_sort) {
    /* Everything is ordered. Just need the right alignment .... */
    output_arena_obj = _align_arena(input_arena_obj, alignment,
				    &start_padding, &end_padding);
    if (!output_arena_obj) {
      return NULL;
    }

    /* ... and to set the popcount indicies */
    set_popcount_indicies(num_fingerprints, num_bits, ordering, popcount_indices);
    
    /* Everything is aligned and ordered, so we're done */
    return Py_BuildValue("iiO", start_padding, end_padding, output_arena_obj);
  }

  /* Not ordered. Make space for the results. */
  output_arena_obj = _alloc_aligned_arena(input_arena_size, alignment,
					  &start_padding, &end_padding);
  if (!output_arena_obj) {
    return NULL;
  }
  output_arena = (unsigned char *)(PyString_AS_STRING(output_arena_obj) + start_padding);

  Py_BEGIN_ALLOW_THREADS;
  qsort(ordering, num_fingerprints, sizeof(ChemFPOrderedPopcount), compare_by_popcount);


  /* Build the new arena based on the values in the old arena */
  for (i=0; i<num_fingerprints; i++) {
    memcpy(output_arena+(i*storage_size), input_arena+(ordering[i].index * storage_size),
	   storage_size);
  }

  /* Create the popcount indicies */
  set_popcount_indicies(num_fingerprints, num_bits, ordering, popcount_indices);


  Py_END_ALLOW_THREADS;

  return Py_BuildValue("iiO", start_padding, end_padding, output_arena_obj);
}


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
  int *popcount_indices, popcount_indices_size;

  if (!PyArg_ParseTuple(args, "iit#iiw#w#",
                        &num_bits,
                        &storage_size, &arena, &arena_size,
                        &start, &end,
                        &ordering, &ordering_size,
                        &popcount_indices, &popcount_indices_size
                        )) {
    return NULL;
  }

  if (bad_num_bits(num_bits) ||
      bad_arena_limits("", arena_size, storage_size, &start, &end) ||
      bad_popcount_indices("", 0, num_bits, popcount_indices_size, NULL)) {
    return NULL;
  }
  if ((ordering_size / sizeof(ChemFPOrderedPopcount)) < (end-start)) {
    PyErr_SetString(PyExc_ValueError, "allocated ordering space is too small");
    return NULL;
  }
  /* TODO: compute the counts first. If everything is in order then */
  /* there's no need to allocate a new arena string. */
  if (end <= start) {
    py_arena = PyString_FromStringAndSize("", 0);
  } else {
    py_arena = PyString_FromStringAndSize(NULL, (end-start)*storage_size);
  }
  if (!py_arena)
    goto error;

  Py_BEGIN_ALLOW_THREADS;
  chemfp_reorder_by_popcount(num_bits, storage_size,
                             arena, start, end,
                             (unsigned char *) PyString_AS_STRING(py_arena),
                             ordering, popcount_indices);
  Py_END_ALLOW_THREADS;

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
  const unsigned char *query_arena, *target_arena;
  int query_start_padding, query_end_padding;
  int query_storage_size, query_arena_size=0, query_start=0, query_end=0;
  int target_start_padding, target_end_padding;
  int target_storage_size, target_arena_size=0, target_start=0, target_end=0;
  
  int *target_popcount_indices, target_popcount_indices_size;
  
  int result_counts_size, *result_counts;

  if (!PyArg_ParseTuple(args, "diiiis#iiiiis#iis#w#",
                        &threshold,
                        &num_bits,
			&query_start_padding, &query_end_padding,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
			&target_start_padding, &target_end_padding,
                        &target_storage_size, &target_arena, &target_arena_size,
                        &target_start, &target_end,
                        &target_popcount_indices, &target_popcount_indices_size,
                        &result_counts, &result_counts_size))
    return NULL;

  if (bad_threshold(threshold) ||
      bad_num_bits(num_bits) ||
      bad_padding("query ", query_start_padding, query_end_padding,
		  &query_arena, &query_arena_size) ||
      bad_padding("target ", target_start_padding, target_end_padding,
		  &target_arena, &target_arena_size) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indices("target ", 1, num_bits, 
                            target_popcount_indices_size, &target_popcount_indices)) {
    return NULL;
  }

  if (query_start > query_end) {
    Py_RETURN_NONE;
  }

  if (result_counts_size < (query_end - query_start)*sizeof(int)) {
    PyErr_SetString(PyExc_ValueError, "not enough space allocated for result_counts");
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS;
  chemfp_count_tanimoto_arena(threshold,
                              num_bits,
                              query_storage_size, query_arena, query_start, query_end,
                              target_storage_size, target_arena, target_start, target_end,
                              target_popcount_indices,
                              result_counts);
  Py_END_ALLOW_THREADS;

  Py_RETURN_NONE;
}
    

/* threshold_tanimoto_arena */
static PyObject *
threshold_tanimoto_arena(PyObject *self, PyObject *args) {
  double threshold;
  int num_bits;
  int query_start_padding, query_end_padding;
  int query_storage_size, query_arena_size, query_start, query_end;
  const unsigned char *query_arena;
  int target_start_padding, target_end_padding;
  int target_storage_size, target_arena_size, target_start, target_end;
  const unsigned char *target_arena;

  int *target_popcount_indices, target_popcount_indices_size;
  int num_cells;
  int *result_offsets, result_offsets_size, result_offsets_start;
  int *result_indices, result_indices_size, result_scores_size;
  double *result_scores;

  int result;

    
  if (!PyArg_ParseTuple(args, "diiiit#iiiiit#iit#w#iw#w#",
                        &threshold,
                        &num_bits,
			&query_start_padding, &query_end_padding,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
			&target_start_padding, &target_end_padding,
                        &target_storage_size, &target_arena, &target_arena_size,
                        &target_start, &target_end,
                        &target_popcount_indices, &target_popcount_indices_size,
                        &result_offsets, &result_offsets_size, &result_offsets_start,
                        &result_indices, &result_indices_size,
                        &result_scores, &result_scores_size
                        ))
    return NULL;

  if (bad_threshold(threshold) ||
      bad_num_bits(num_bits) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_padding("query ", query_start_padding, query_end_padding, 
		  &query_arena, &query_arena_size) ||
      bad_padding("target ", target_start_padding, target_end_padding, 
		  &target_arena, &target_arena_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indices("target ", 1, num_bits,
                            target_popcount_indices_size, &target_popcount_indices)) {
    return NULL;
  }

  if (bad_offsets(query_end-query_start, result_offsets_size, result_offsets_start) ||
      bad_cells(target_end-target_start, result_indices_size,
                result_scores_size, &num_cells)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS;
  result = chemfp_threshold_tanimoto_arena(
        threshold,
        num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        target_popcount_indices,
        result_offsets+result_offsets_start,
        num_cells, result_indices, result_scores);
  Py_END_ALLOW_THREADS;

  return PyInt_FromLong(result);
}

/* knearest_tanimoto_arena */
static PyObject *
knearest_tanimoto_arena(PyObject *self, PyObject *args) {
  int k;
  double threshold;
  int num_bits;
  int query_start_padding, query_end_padding;
  int query_storage_size, query_arena_size, query_start, query_end;
  const unsigned char *query_arena;
  int target_start_padding, target_end_padding;
  int target_storage_size, target_arena_size, target_start, target_end;
  const unsigned char *target_arena;

  int *target_popcount_indices, target_popcount_indices_size;
  int num_cells;
  int *result_offsets, result_offsets_size, result_offsets_start;
  int *result_indices, result_indices_size, result_scores_size;
  double *result_scores;

  int result;

    
  if (!PyArg_ParseTuple(args, "idiiiit#iiiiit#iit#w#iw#w#",
                        &k, &threshold,
                        &num_bits,
			&query_start_padding, &query_end_padding,
                        &query_storage_size, &query_arena, &query_arena_size,
                        &query_start, &query_end,
			&target_start_padding, &target_end_padding,
                        &target_storage_size, &target_arena, &target_arena_size,
                        &target_start, &target_end,
                        &target_popcount_indices, &target_popcount_indices_size,
                        &result_offsets, &result_offsets_size, &result_offsets_start,
                        &result_indices, &result_indices_size,
                        &result_scores, &result_scores_size
                        ))
    return NULL;

  if (bad_k(k) ||
      bad_threshold(threshold) ||
      bad_num_bits(num_bits) ||
      bad_padding("query ", query_start_padding, query_end_padding,
		  &query_arena, &query_arena_size) ||
      bad_padding("target ", target_start_padding, target_end_padding,
		  &target_arena, &target_arena_size) ||
      bad_fingerprint_sizes(num_bits, query_storage_size, target_storage_size) ||
      bad_arena_limits("query ", query_arena_size, query_storage_size,
		       &query_start, &query_end) ||
      bad_arena_limits("target ", target_arena_size, target_storage_size,
		       &target_start, &target_end) ||
      bad_popcount_indices("target ", 1, num_bits,
                            target_popcount_indices_size, &target_popcount_indices)) {
    return NULL;
  }

  if (bad_offsets(query_end-query_start, result_offsets_size, result_offsets_start) ||
      bad_cells(k, result_indices_size, result_scores_size, &num_cells)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS;
  result = chemfp_knearest_tanimoto_arena(
        k, threshold,
        num_bits,
        query_storage_size, query_arena, query_start, query_end,
        target_storage_size, target_arena, target_start, target_end,
        target_popcount_indices,
        result_offsets+result_offsets_start,
        num_cells, result_indices, result_scores);
  Py_END_ALLOW_THREADS;

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

  /* FPS */
  {"fps_line_validate", fps_line_validate, METH_VARARGS,
   "fps_line_validate (TODO: document)"},
  {"fps_parse_id_fp", fps_parse_id_fp, METH_VARARGS,
   "fps_parse_id_fp (TODO: document)"},

  {"fps_threshold_tanimoto_search", fps_threshold_tanimoto_search, METH_VARARGS,
   "fps_threshold_tanimoto_search (TODO: document)"},

  {"fps_count_tanimoto_hits", fps_count_tanimoto_hits, METH_VARARGS,
   "fps_count_tanimoto_hits (TODO: document)"},


  {"fps_knearest_search_init", fps_knearest_search_init, METH_VARARGS,
   "fps_knearest_search_init (TODO: document)"},
  {"fps_knearest_tanimoto_search_feed", fps_knearest_tanimoto_search_feed, METH_VARARGS,
   "fps_knearest_tanimoto_search_feed (TODO: document)"},
  {"fps_knearest_search_finish", fps_knearest_search_finish, METH_VARARGS,
   "fps_knearest_search_finish (TODO: document)"},
  {"fps_knearest_search_free", fps_knearest_search_free, METH_VARARGS,
   "fps_knearest_search_free (TODO: document)"},

  {"count_tanimoto_arena", count_tanimoto_arena, METH_VARARGS,
   "count_tanimoto_arena (TODO: document)"},

  {"threshold_tanimoto_arena", threshold_tanimoto_arena, METH_VARARGS,
   "threshold_tanimoto_arena (TODO: document)"},

  {"knearest_tanimoto_arena", knearest_tanimoto_arena, METH_VARARGS,
   "knearest_tanimoto_arena (TODO: document)"},

  {"reorder_by_popcount", reorder_by_popcount, METH_VARARGS,
   "reorder_by_popcount (TODO: document)"},

  {"make_sorted_aligned_arena", make_sorted_aligned_arena, METH_VARARGS,
   "make_sorted_aligned_arena (TODO: document)"},
  {"make_unsorted_aligned_arena", make_unsorted_aligned_arena, METH_VARARGS,
   "make_unsorted_aligned_arena (TODO: document)"},

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
