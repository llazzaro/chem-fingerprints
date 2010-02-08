#include <Python.h>

#include "chemfp.h"

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

static PyMethodDef chemfp_methods[] = {
  {"hex_isvalid", hex_isvalid, METH_VARARGS, "is this a valid hex fingerprint"},
  {"hex_popcount", hex_popcount, METH_VARARGS, "popcount"},
  {"hex_intersect_popcount", hex_intersect_popcount, METH_VARARGS, "intersect_popcount"},
  {"hex_tanimoto", hex_tanimoto, METH_VARARGS, "Tanimoto"},
  {"hex_contains", hex_contains, METH_VARARGS, "contains"},


  {"byte_popcount", byte_popcount, METH_VARARGS, "popcount"},
  {"byte_intersect_popcount", byte_intersect_popcount, METH_VARARGS, "intersect_popcount"},
  {"byte_tanimoto", byte_tanimoto, METH_VARARGS, "Tanimoto"},
  {"byte_contains", byte_contains, METH_VARARGS, "contains"},

  {"nlargest_tanimoto_block", nlargest_tanimoto_block, METH_VARARGS,
   "nlargest_tanimoto_block"},

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
