#include <float.h>

#include "pysearch_results.h"
#include "structmember.h"

#include "chemfp_internal.h"

/************ Search Result type ***************/


/* Help with cyclical garbage collection, in case someone does result.target_ids = result */
static int
SearchResults_traverse(SearchResults *self, visitproc visit, void *arg) {
  Py_VISIT(self->target_ids);
  return 0;
}

static int
SearchResults_clear_memory(SearchResults *self) {
  if (self->results) {
    chemfp_free_results(self->num_results, self->results);
    self->results = NULL;
  }
  self->num_results = 0;

  Py_CLEAR(self->target_ids);
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
    self->target_ids = Py_None;
    return (PyObject *)self;
}

static int
SearchResults_init(SearchResults *self, PyObject *args, PyObject *kwds)
{
  int num_results=0;
  PyObject *target_ids=Py_None;
  chemfp_search_result *results;

  static char *kwlist[] = {"num_results", "target_ids", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|O", kwlist, &num_results, &target_ids)) {
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

  Py_XINCREF(target_ids);
  Py_XDECREF(self->target_ids);
  self->target_ids = target_ids;

  return 0;
}

/* len(search_results) */
static Py_ssize_t
SearchResults_length(SearchResults *result) {
  return (Py_ssize_t) result->num_results;
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

static int
check_min_max_score(PyObject *min_score_obj, PyObject *max_score_obj,
                    double *min_score, double *max_score) {
  double value;
  if (min_score_obj == Py_None) {
    *min_score = -HUGE_VAL;
  } else {
    value = PyFloat_AsDouble(min_score_obj);
    if (value == -1.0) {
      if (PyErr_Occurred()) {
        return 0;
      }
    }
    *min_score = value;
  }
  if (max_score_obj == Py_None) {
    *max_score = HUGE_VAL;
  } else {
    value = PyFloat_AsDouble(max_score_obj);
    if (value == -1.0) {
      if (PyErr_Occurred()) {
        return 0;
      }
    }
    *max_score = value;
  }
  return 1;
}

static int
check_interval(const char *interval, int *include_min, int *include_max) {
  switch (interval[0]) {
  case '(':
    *include_min = 0;
    break;
  case '[':
    *include_min = 1;
    break;
  default:
      PyErr_SetString(PyExc_ValueError, "First interval character must be '(' or '['");
      return 0;
  }
  switch (interval[1]) {
  case ')':
    *include_max = 0;
    break;
  case ']':
    *include_max = 1;
    break;
  default:
      PyErr_SetString(PyExc_ValueError, "Second interval character must be ')' or ']'");
      return 0;
  }
  if (interval[2]) {
    PyErr_SetString(PyExc_ValueError, "The interval may only contain two characters");
    return 0;
  }
  return 1;
}


static PyObject *
SearchResults_reorder_all(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"order", NULL};
  const char *ordering = "decreasing-score";
  int errval;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s:reorder_all", kwlist, &ordering)) {
    return NULL;
  }
  errval = chemfp_search_results_reorder(self->num_results, self->results, ordering);
  if (errval) {
    PyErr_SetString(PyExc_ValueError, chemfp_strerror(errval));
    return NULL;
  }
  Py_RETURN_NONE;
}

static PyObject *
SearchResults_reorder_row(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", "order", NULL};
  int row=-1;
  int errval;
  const char *ordering = "decreasing-score";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|s:reorder_row", kwlist, &row, &ordering)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  errval = chemfp_search_result_reorder(self->results+row, ordering);
  if (errval) {
    PyErr_SetString(PyExc_ValueError, chemfp_strerror(errval));
    return NULL;
  }
  Py_RETURN_NONE;
}

#define COUNT_ALL_MACRO(expr)                           \
  for (row=0; row<num_rows; row++) {                    \
    num_hits = chemfp_get_num_hits(self->results+row);  \
    scores = self->results[row].scores;                 \
    for (i=0; i<num_hits; i++) {                        \
      if (expr) {                                       \
        count++;                                        \
      }                                                 \
    }                                                   \
  }

static PyObject *
SearchResults_count_all(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"min_score", "max_score", "interval", NULL};
  int row=-1;
  PyObject *min_score_obj=Py_None, *max_score_obj=Py_None;
  double min_score=0.0, max_score=1.0;
  int num_hits, num_rows;
  double *scores;
  int include_min=0, include_max=0;
  int count=0;
  int i;
  const char *interval = "[]";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOs:count_all", kwlist,
                                   &min_score_obj, &max_score_obj, &interval)) {
    return NULL;
  }
  if (!check_min_max_score(min_score_obj, max_score_obj, &min_score, &max_score) ||
      !check_interval(interval, &include_min, &include_max)) {
    return NULL;
  }
  if ((min_score > max_score) ||
      (min_score == max_score && !(include_min && include_max))) {
    return PyInt_FromLong(0);
  }

  num_rows = self->num_results;
  if (include_min) {
    if (include_max) {
      /* [] -- Include both ends */
      if (min_score_obj == Py_None) {
        if (max_score_obj == Py_None) {
          /* Special case; just report the number of elements */
          for (row=0; row<num_rows; row++) {
            count += chemfp_get_num_hits(self->results+row);
          }
        } else {
          /* No lower bound */
          COUNT_ALL_MACRO(scores[i] <= max_score);
        }
      } else {
        if (max_score_obj == Py_None) {
          /* Lower bound but no upper bound */
          COUNT_ALL_MACRO(min_score <= scores[i]);
        } else {
          /* Definite lower and upper bound */
          COUNT_ALL_MACRO(min_score <= scores[i] && scores[i] <= max_score);
        }
      }
    } else {
      /* [) -- Include the minimum but not the maximum */
      if (min_score_obj == Py_None) {
        /* There is no minimum */
        COUNT_ALL_MACRO(scores[i] < max_score);
      } else {
        /* There is a minimum and a maximum test */
        COUNT_ALL_MACRO(min_score <= scores[i] && scores[i] < max_score);
      }
    }
  } else {
    if (include_max) {
      /* (] -- Exclude the minimum, include the maximum */
      if (max_score_obj == Py_None) {
        /* No specified max, so must only be greater than the lower bound */
        COUNT_ALL_MACRO(min_score < scores[i]);
      } else {
        COUNT_ALL_MACRO(min_score < scores[i] && scores[i] <= max_score);
      }
    } else {
      /* () -- Exclude the minimum and exclude the maximum */
      COUNT_ALL_MACRO(min_score < scores[i] && scores[i] < max_score);
    }
  }
  return PyInt_FromLong(count);
}

#define COUNT_ROW_MACRO(expr)                         \
  for (i=0; i<num_hits; i++) {                        \
    if (expr) {                                       \
      count++;                                        \
    }                                                 \
  }

static PyObject *
SearchResults_count_row(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", "min_score", "max_score", "interval", NULL};
  int row=-1;
  PyObject *min_score_obj=Py_None, *max_score_obj=Py_None;
  double min_score=0.0, max_score=1.0;
  int num_hits;
  double *scores;
  int include_min=0, include_max=0;
  int count=0;
  int i;
  const char *interval = "[]";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|OOs:count_row", kwlist,
                                   &row, &min_score_obj, &max_score_obj, &interval)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row) ||
      !check_min_max_score(min_score_obj, max_score_obj, &min_score, &max_score) ||
      !check_interval(interval, &include_min, &include_max)) {
    return NULL;
  }
  num_hits = chemfp_get_num_hits(self->results+row);
  scores = (self->results+row)->scores;

  if ((min_score > max_score) ||
      (min_score == max_score && (!include_min || !include_max))) {
    return PyInt_FromLong(0);
  }

  if (include_min) {
    if (include_max) {
      /* [] -- Include both ends */
      if (min_score_obj == Py_None) {
        if (max_score_obj == Py_None) {
          /* Special case; just report the number of elements */
          count = chemfp_get_num_hits(self->results+row);
        } else {
          /* No lower bound */
          COUNT_ROW_MACRO(scores[i] <= max_score);
        }
      } else {
        if (max_score_obj == Py_None) {
          /* Lower bound but no upper bound */
          COUNT_ROW_MACRO(min_score <= scores[i]);
        } else {
          /* Definite lower and upper bound */
          COUNT_ROW_MACRO(min_score <= scores[i] && scores[i] <= max_score);
        }
      }
    } else {
      /* [) -- Include the minimum but not the maximum */
      if (min_score_obj == Py_None) {
        /* There is no minimum */
        COUNT_ROW_MACRO(scores[i] < max_score);
      } else {
        /* There is a minimum and a maximum test */
        COUNT_ROW_MACRO(min_score <= scores[i] && scores[i] < max_score);
      }
    }
  } else {
    if (include_max) {
      /* (] -- Exclude the minimum, include the maximum */
      if (max_score_obj == Py_None) {
        /* No specified max, so must only be greater than the lower bound */
        COUNT_ROW_MACRO(min_score < scores[i]);
      } else {
        COUNT_ROW_MACRO(min_score < scores[i] && scores[i] <= max_score);
      }
    } else {
      /* () -- Exclude the minimum and exclude the maximum */
      COUNT_ROW_MACRO(min_score < scores[i] && scores[i] < max_score);
    }
  }
  return PyInt_FromLong(count);
}

#define CUMULATIVE_SCORE_ALL_MACRO(expr)                \
  for (row=0; row<num_rows; row++) {                    \
    num_hits = chemfp_get_num_hits(self->results+row);  \
    scores = self->results[row].scores;                 \
    for (i=0; i<num_hits; i++) {                        \
      if (expr) {                                       \
        score += scores[i];                             \
      }                                                 \
    }                                                   \
  }

static PyObject *
SearchResults_cumulative_score_all(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"min_score", "max_score", "interval", NULL};
  int row=-1;
  PyObject *min_score_obj=Py_None, *max_score_obj=Py_None;
  double min_score=0.0, max_score=1.0;
  int num_hits, num_rows;
  double *scores;
  int include_min=0, include_max=0;
  double score=0.0;
  int i;
  const char *interval = "[]";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOs:cumulative_score_all", kwlist,
                                   &min_score_obj, &max_score_obj, &interval)) {
    return NULL;
  }
  if (!check_min_max_score(min_score_obj, max_score_obj, &min_score, &max_score) ||
      !check_interval(interval, &include_min, &include_max)) {
    return NULL;
  }
  if ((min_score > max_score) ||
      (min_score == max_score && !(include_min && include_max))) {
    return PyInt_FromLong(0);
  }

  num_rows = self->num_results;
  if (include_min) {
    if (include_max) {
      /* [] -- Include both ends */
      if (min_score_obj == Py_None) {
        if (max_score_obj == Py_None) {
          CUMULATIVE_SCORE_ALL_MACRO(1);
        } else {
          /* No lower bound */
          CUMULATIVE_SCORE_ALL_MACRO(scores[i] <= max_score);
        }
      } else {
        if (max_score_obj == Py_None) {
          /* Lower bound but no upper bound */
          CUMULATIVE_SCORE_ALL_MACRO(min_score <= scores[i]);
        } else {
          /* Definite lower and upper bound */
          CUMULATIVE_SCORE_ALL_MACRO(min_score <= scores[i] && scores[i] <= max_score);
        }
      }
    } else {
      /* [) -- Include the minimum but not the maximum */
      if (min_score_obj == Py_None) {
        /* There is no minimum */
        CUMULATIVE_SCORE_ALL_MACRO(scores[i] < max_score);
      } else {
        /* There is a minimum and a maximum test */
        CUMULATIVE_SCORE_ALL_MACRO(min_score <= scores[i] && scores[i] < max_score);
      }
    }
  } else {
    if (include_max) {
      /* (] -- Exclude the minimum, include the maximum */
      if (max_score_obj == Py_None) {
        /* No specified max, so must only be greater than the lower bound */
        CUMULATIVE_SCORE_ALL_MACRO(min_score < scores[i]);
      } else {
        CUMULATIVE_SCORE_ALL_MACRO(min_score < scores[i] && scores[i] <= max_score);
      }
    } else {
      /* () -- Exclude the minimum and exclude the maximum */
      CUMULATIVE_SCORE_ALL_MACRO(min_score < scores[i] && scores[i] < max_score);
    }
  }
  return PyFloat_FromDouble(score);
}

#define CUMULATIVE_SCORE_ROW_MACRO(expr)              \
  for (i=0; i<num_hits; i++) {                        \
    if (expr) {                                       \
      score += scores[i];                             \
    }                                                 \
  }


static PyObject *
SearchResults_cumulative_score_row(SearchResults *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"row", "min_score", "max_score", "interval", NULL};
  int row=-1;
  PyObject *min_score_obj=Py_None, *max_score_obj=Py_None;
  double min_score=0.0, max_score=1.0;
  int num_hits;
  double *scores;
  int include_min=0, include_max=0;
  double score=0.0;
  int i;
  const char *interval = "[]";
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|OOs:cumulative_score_row", kwlist,
                                   &row, &min_score_obj, &max_score_obj, &interval)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row) ||
      !check_min_max_score(min_score_obj, max_score_obj, &min_score, &max_score) ||
      !check_interval(interval, &include_min, &include_max)) {
    return NULL;
  }
  num_hits = chemfp_get_num_hits(self->results+row);
  scores = (self->results+row)->scores;

  if ((min_score > max_score) ||
      (min_score == max_score && (!include_min || !include_max))) {
    return PyFloat_FromDouble(0.0);
  }

  if (include_min) {
    if (include_max) {
      /* [] -- Include both ends */
      if (min_score_obj == Py_None) {
        if (max_score_obj == Py_None) {
          CUMULATIVE_SCORE_ROW_MACRO(1);
        } else {
          /* No lower bound */
          CUMULATIVE_SCORE_ROW_MACRO(scores[i] <= max_score);
        }
      } else {
        if (max_score_obj == Py_None) {
          /* Lower bound but no upper bound */
          CUMULATIVE_SCORE_ROW_MACRO(min_score <= scores[i]);
        } else {
          /* Definite lower and upper bound */
          CUMULATIVE_SCORE_ROW_MACRO(min_score <= scores[i] && scores[i] <= max_score);
        }
      }
    } else {
      /* [) -- Include the minimum but not the maximum */
      if (min_score_obj == Py_None) {
        /* There is no minimum */
        CUMULATIVE_SCORE_ROW_MACRO(scores[i] < max_score);
      } else {
        /* There is a minimum and a maximum test */
        CUMULATIVE_SCORE_ROW_MACRO(min_score <= scores[i] && scores[i] < max_score);
      }
    }
  } else {
    if (include_max) {
      /* (] -- Exclude the minimum, include the maximum */
      if (max_score_obj == Py_None) {
        /* No specified max, so must only be greater than the lower bound */
        CUMULATIVE_SCORE_ROW_MACRO(min_score < scores[i]);
      } else {
        CUMULATIVE_SCORE_ROW_MACRO(min_score < scores[i] && scores[i] <= max_score);
      }
    } else {
      /* () -- Exclude the minimum and exclude the maximum */
      CUMULATIVE_SCORE_ROW_MACRO(min_score < scores[i] && scores[i] < max_score);
    }
  }
  return PyFloat_FromDouble(score);
}


static PyObject *
SearchResults_clear_all(SearchResults *self) {
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
SearchResults_get_indices_and_scores(SearchResults *self, PyObject *args, PyObject *kwds) {
  int n, i;
  PyObject *hits, *obj;
  chemfp_search_result *result;
  static char *kwlist[] = {"row", NULL};
  int row;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:get_indices", kwlist, &row)) {
    return NULL;
  }
  if (!check_row(self->num_results, &row)) {
    return NULL;
  }
  result = self->results+row;
  n = chemfp_get_num_hits(result);
  
  hits = PyList_New(n);
  if (!hits) {
    return NULL;
  }
  for (i=0; i<n; i++) {
    obj = Py_BuildValue("(id)", result->indices[i], result->scores[i]);
    if (!obj) {
      goto error;
    }
    PyList_SET_ITEM(hits, i, obj);
  }
  return hits;

 error:
  /* I could be smarter and only deallocate the ones which I know are present. */
  /* Proving that's correct is a bit trickier than I want to do. */
  for (i=0; i<n; i++) {
    Py_XDECREF(PyList_GET_ITEM(hits, i));
  }
  Py_DECREF(hits);
  return NULL;
}


static PyObject *
new_array(const char *datatype) {
  PyObject *array_module;
  static PyObject *array_new = NULL;
  if (!array_new) {
    array_module = PyImport_ImportModule("array");
    if (!array_module) {
      return NULL;
    }
    array_new = PyObject_GetAttrString(array_module, "array");
    if (!array_new) {
      return NULL;
    }
  }
  return PyObject_CallFunction(array_new, "s", datatype);
}

static PyObject *
data_blob_to_array(int num_hits, void *data_blob, const char *array_type, size_t num_bytes_per_item) {
  PyObject *array=NULL, *ignore=NULL, *fromstring=NULL;
  
  if (!num_hits) {
    /* Special case: return a tuple if there are no hits */
    /* (I'm not sure if this is a good idea, but it's fast. I justify it */
    /* by saying that this function returns back something which supports */
    /* len(), [index], and forward/reverse iteration */
    return PyTuple_New(0);
  }
  array = new_array(array_type);
  if (!array) {
    goto error;
  }
  fromstring = PyObject_GetAttrString(array, "fromstring");
  if (!fromstring) {
    goto error;
  }
  ignore = PyObject_CallFunction(fromstring, "s#", data_blob, num_hits*num_bytes_per_item);
  if (!ignore) {
    goto error;
  }

  Py_DECREF(ignore);
  Py_DECREF(fromstring);
  return array;

 error:
  if (fromstring) {
    Py_DECREF(fromstring);
  }
  if (array) {
    Py_DECREF(array);
  }
  return NULL;
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
  return data_blob_to_array(chemfp_get_num_hits(self->results+row),
                            self->results[row].indices,
                            "i", sizeof(int));
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
  return data_blob_to_array(chemfp_get_num_hits(self->results+row),
                            self->results[row].scores,
                            "d", sizeof(double));
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
  {"clear_all", (PyCFunction) SearchResults_clear_all, METH_VARARGS | METH_KEYWORDS,
   "Removes all hits from all of the search results"},
  {"_clear_row", (PyCFunction) SearchResults_clear_row, METH_VARARGS | METH_KEYWORDS,
   "(internal) Remove all hits from a given row result"},
  {"cumulative_score_all", (PyCFunction) SearchResults_cumulative_score_all, METH_VARARGS | METH_KEYWORDS,
   "The sum of all scores in all rows which are between `min_score` and `max_score`"},
  {"_cumulative_score_row", (PyCFunction) SearchResults_cumulative_score_row, METH_VARARGS | METH_KEYWORDS,
   "(internal) The sum of the scores which are between `min_score` and `max_score` for a given row"},
  {"count_all", (PyCFunction) SearchResults_count_all, METH_VARARGS | METH_KEYWORDS,
   "Count the number of hits with a score between `min_score` and `max_score`"},
  {"_count_row", (PyCFunction) SearchResults_count_row, METH_VARARGS | METH_KEYWORDS,
   "(internal) Count the number of hits with a score between `min_score` and `max_score` for a given row"},
  {"_get_indices", (PyCFunction) SearchResults_get_indices, METH_VARARGS | METH_KEYWORDS,
   "(internal) The list of target indices, in the current ordering for a given row"},
  {"_get_scores", (PyCFunction) SearchResults_get_scores, METH_VARARGS | METH_KEYWORDS,
   "(internal) The list of target scores, in the current ordering for a given row"},
  {"_get_indices_and_scores", (PyCFunction) SearchResults_get_indices_and_scores, METH_VARARGS | METH_KEYWORDS,
   "(internal) The list of (target index, score) pairs, in the current ordering for a given row"},
  {"_size", (PyCFunction) SearchResults_size, METH_VARARGS | METH_KEYWORDS,
   "(internal) The number of hits for a given row"},
  {"reorder_all", (PyCFunction) SearchResults_reorder_all, METH_VARARGS | METH_KEYWORDS,
   "Reorder the hits for all of the rows based on the requested ordering"},
  {"_reorder_row", (PyCFunction) SearchResults_reorder_row, METH_VARARGS | METH_KEYWORDS,
   "(internal) Reorder the hits based on the requested ordering for a given row"},
  {"_add_hit", (PyCFunction) SearchResults_add_hit, METH_VARARGS | METH_KEYWORDS,
   "(internal) Add a target index and hit score to a given row"},
  {NULL}
};

static PyMemberDef SearchResults_members[] = {
  {"target_ids", T_OBJECT_EX, offsetof(SearchResults, target_ids), 0,
   "list of fingerprint identifiers"},
  {NULL}
};


static PySequenceMethods SearchResults_as_sequence = {
    (lenfunc)SearchResults_length,                       /* sq_length */
    NULL,       /* sq_concat */
    NULL,       /* sq_repeat */
    NULL,       /* sq_item */
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
