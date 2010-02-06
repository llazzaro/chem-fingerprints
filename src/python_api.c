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
nlargest_tanimoto(PyObject *self, PyObject *args) {
  unsigned char *query_fp, *fps, *indicies_buffer, *scores_buffer;
  int query_len, fps_len, num_fps, offset, storage_len;
  int indicies_len, scores_len;
  int n;
  if (!PyArg_ParseTuple(args, "s#s#iiiw#w#",
						&query_fp, &query_len,
						&fps, &fps_len,
						&offset, &storage_len,
						&n, &indicies_buffer, &indicies_len,
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

  if (offset > fps_len) {
	PyErr_SetString(PyExc_TypeError,
					"offset is larger than the fps buffer");
	return NULL;
  }
  fps_len -= offset;
  fps += offset;
  offset = 0;
  if (fps_len == 0)
	return PyInt_FromLong(0);

  if ((fps_len % storage_len) != 0) {
	PyErr_SetString(PyExc_TypeError,
					"fps length is not a multiple of the storage size");
	return NULL;
  }
  num_fps = fps_len / storage_len;
  if (n > num_fps)
	n = num_fps;

  if (n * sizeof(int) > indicies_len) {
	PyErr_SetString(PyExc_TypeError, "indicies buffer is not long enough");
	return NULL;
  }
  if (n * sizeof(double) > scores_len) {
	PyErr_SetString(PyExc_TypeError, "score buffer is not long enough");
	return NULL;
  }

  return PyInt_FromLong(chemfp_nlargest_tanimoto(
		query_len, query_fp,		
		num_fps, fps, offset, storage_len,
		n, (int *) indicies_buffer, (double *) scores_buffer));
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

  {"nlargest_tanimoto", nlargest_tanimoto, METH_VARARGS, "nlargest_tanimoto"},

  {NULL, NULL, 0, NULL}        /* Sentinel */

};

PyMODINIT_FUNC
init_chemfp(void)
{
    (void) Py_InitModule("_chemfp", chemfp_methods);
}
