#include "pysearch_results.h"
#include "structmember.h"

#include "chemfp_internal.h"

/************ Search Result type ***************/


/* Help with cyclical garbage collection, in case someone does result.ids = result */
static int
SearchResults_traverse(SearchResults *self, visitproc visit, void *arg) {
  Py_VISIT(self->ids);
  return 0;
}

static int
SearchResults_clear_memory(SearchResults *self) {
  if (self->results) {
    chemfp_free_results(self->num_results, self->results);
    self->results = NULL;
  }
  self->num_results = 0;

  Py_CLEAR(self->ids);
  return 0;
}

static void
SearchResults_dealloc(SearchResults *self) {
  SearchResults_clear_memory(self);
  self->ob_type->tp_free((PyObject *) self);
}

static PyObject *
SearchResults_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    SearchResults *self;

    self = (SearchResults *) type->tp_alloc(type, 0);
    if (self == NULL) {
      return NULL;
    }
    self->num_results = 0;
    self->results = NULL;
    Py_INCREF(Py_None);
    self->ids = Py_None;
    return (PyObject *)self;
}

static int
SearchResults_init(SearchResults *self, PyObject *args, PyObject *kwds)
{
  int num_results=0;
  PyObject *ids=Py_None;
  chemfp_search_result *results;

  static char *kwlist[] = {"num_results", "ids", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|O", kwlist, &num_results, &ids)) {
    return -1;
  }
  if (num_results < 0) {
    PyErr_SetString(PyExc_ValueError, "num_results must be non-negative");
    return -1;
  }
  if (num_results == 0) {
    results = NULL;
  } else {
    results = chemfp_alloc_search_results(num_results);
    if (!results) {
      PyErr_NoMemory();
      return -1;
    }
  }

  self->num_results = num_results;
  self->results = results;

  Py_XINCREF(ids);
  Py_XDECREF(self->ids);
  self->ids = ids;

  return 0;
}

static Py_ssize_t
SearchResults_length(SearchResults *result) {
  return (Py_ssize_t) result->num_results;
}

static int assign_threshold_hits(void *data, int i, int target_index, double score) {
  PyObject *list = (PyObject *) data;
  PyObject *tuple;
  
  tuple = Py_BuildValue("(id)", target_index, score);
  if (!tuple) {
    return 1;
  }
  PyList_SET_ITEM(list, i, tuple);
  return 0;
}

static int assign_threshold_indices(void *data, int i, int target_index, double score) {
  PyObject *list = (PyObject *) data;
  PyObject *index;
  
  index = PyInt_FromLong(target_index);
  if (!index) {
    return 1;
  }
  PyList_SET_ITEM(list, i, index);
  return 0;
}

static int assign_threshold_scores(void *data, int i, int target_index, double score) {
  PyObject *list = (PyObject *) data;
  PyObject *pyscore;
  
  pyscore = PyFloat_FromDouble(score);
  if (!pyscore) {
    return 1;
  }
  PyList_SET_ITEM(list, i, pyscore);
  return 0;
}

static PyObject *
_fill_row(chemfp_search_result *result, chemfp_assign_hits_p extract_func) {
  int errval, n, i;
  PyObject *hits, *obj;

  n = chemfp_get_num_hits(result);
  hits = PyList_New(n);
  if (!hits) {
    return NULL;
  }
  errval = chemfp_search_result_get_hits(result, extract_func, hits);
  if (errval) {
    for (i=0; i<n; i++) {
      obj = PyList_GetItem(hits, i);
      if (obj) {
        Py_DECREF(obj);
      } else {
        break;
      }
    }
    Py_DECREF(hits);
    return NULL;
  }
  return hits;
}

static PyObject *
SearchResults_item(SearchResults *self, Py_ssize_t index) {
  if (index < 0 || index >= self->num_results) {
    PyErr_SetString(PyExc_IndexError, "row index is out of range");
    return NULL;
  }
  return _fill_row(self->results+index, assign_threshold_hits);
}


static int
check_row(int num_results, int *row) {
  int row_ = *row;
  if (row_ < 0) {
    row_ = num_results + row_;
    if (row_ < 0) {
      PyErr_SetString(PyExc_IndexError, "row index is out of range");
      return 0;
    }
    *row = row_;
  } else if (row_ >= num_results) {
    PyErr_SetString(PyExc_IndexError, "row index is out of range");
    return 0;
  }
  return 1;
}

static PyObject *
SearchResults_sort(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"order", NULL};
  const char *sort_order = "decreasing";
  int errval;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s:sort", kwlist, &sort_order)) {
    return NULL;
  }
  errval = chemfp_search_results_sort(self->num_results, self->results, sort_order);
  if (errval) {
    PyErr_SetString(PyExc_ValueError, chemfp_strerror(errval));
    return NULL;
  }
  Py_RETURN_NONE;
}

static PyObject *
SearchResults_sort_row(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"order", "row", NULL};
  int row=-1;
  int errval;
  const char *sort_order = "decreasing";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|s:sort", kwlist, &row, &sort_order)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  errval = chemfp_search_result_sort(self->results+row, sort_order);
  if (errval) {
    PyErr_SetString(PyExc_ValueError, chemfp_strerror(errval));
    return NULL;
  }
  Py_RETURN_NONE;
}


static PyObject *
SearchResults_clear(SearchResults *self) {
  int i;
  for (i=0; i<self->num_results; i++) {
    chemfp_search_result_clear(self->results+i);
  }
  Py_RETURN_NONE;
}

static PyObject *
SearchResults_clear_row(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", NULL};
  int row;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:clear", kwlist, &row)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  chemfp_search_result_clear(self->results+row);
  Py_RETURN_NONE;
}


static PyObject *
SearchResults_get_indices(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", NULL};
  int row;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:get_indices", kwlist, &row)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  return _fill_row(self->results+row, assign_threshold_indices);
}

static PyObject *
SearchResults_get_scores(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", NULL};
  int row;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:get_scores", kwlist, &row)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  return _fill_row(self->results+row, assign_threshold_scores);
}

static PyObject *
SearchResults_size(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", NULL};
  int row;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:size", kwlist, &row)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  return PyInt_FromLong(chemfp_get_num_hits(self->results+row));
}


static PyObject *
SearchResults_add_hit(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", "column", "score", NULL};
  int row, column;
  double score;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "iid:_add_hit", kwlist,
                                   &row, &column, &score)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  return PyInt_FromLong(chemfp_add_hit(self->results+row, column, score));
  Py_RETURN_NONE;
}

static PyMethodDef SearchResults_methods[] = {
  {"clear", (PyCFunction) SearchResults_clear, METH_VARARGS | METH_KEYWORDS,
   "clear the hits in-place"},
  {"clear_row", (PyCFunction) SearchResults_clear_row, METH_VARARGS | METH_KEYWORDS,
   "clear the hits in-place"},
  {"get_indices", (PyCFunction) SearchResults_get_indices, METH_VARARGS | METH_KEYWORDS,
   "get the hit indices for a given row"},
  {"get_scores", (PyCFunction) SearchResults_get_scores, METH_VARARGS | METH_KEYWORDS,
   "get the hit scores for a given row"},
  {"size", (PyCFunction) SearchResults_size, METH_VARARGS | METH_KEYWORDS,
   "the number of hits in a given row"},
  {"sort", (PyCFunction) SearchResults_sort, METH_VARARGS | METH_KEYWORDS,
   "Sort search results rows in-place"},
  {"sort_row", (PyCFunction) SearchResults_sort_row, METH_VARARGS | METH_KEYWORDS,
   "Sort search results rows in-place"},
  {"_add_hit", (PyCFunction) SearchResults_add_hit, METH_VARARGS | METH_KEYWORDS,
   "(private method) add a hit"},
  
  {NULL}
};

static PyMemberDef SearchResults_members[] = {
  {"ids", T_OBJECT_EX, offsetof(SearchResults, ids), 0, "list of fingerprint identifiers"},
  {NULL}
};


static PySequenceMethods SearchResults_as_sequence = {
    (lenfunc)SearchResults_length,                       /* sq_length */
    NULL,       /* sq_concat */
    NULL,       /* sq_repeat */
    (ssizeargfunc)SearchResults_item,                    /* sq_item */
    NULL,       /* sq_slice */
    NULL,       /* sq_ass_item */
    NULL,       /* sq_ass_slice */
    NULL,       /* sq_contains */
    NULL,       /* sq_inplace_concat */
    NULL        /* sq_inplace_repeat */
};



PyTypeObject chemfp_py_SearchResultsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "chemfp.search.SearchResults",       /*tp_name*/
    sizeof(SearchResults), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor) SearchResults_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    &SearchResults_as_sequence,          /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    "Documentation goes here",                /* tp_doc */
    (traverseproc) SearchResults_traverse,    /* tp_traverse */
    (inquiry) SearchResults_clear_memory,     /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    SearchResults_methods,     /* tp_methods */
    SearchResults_members,     /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)SearchResults_init,      /* tp_init */
    0,                         /* tp_alloc */
    SearchResults_new          /* tp_new */
};
