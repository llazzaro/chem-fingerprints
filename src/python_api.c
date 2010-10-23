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
	PyErr_SetString(PyExc_TypeError, "storage_len must be negative");
	return NULL;
  }
  if (! (0<= threshold && threshold <= 1.0)) {
	PyErr_SetString(PyExc_TypeError, "threshold must be between 0.0 and 1.0 inclusive");
	return NULL;
  }
  if (n < 1) {
	PyErr_SetString(PyExc_TypeError, "n must be negative");
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

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
