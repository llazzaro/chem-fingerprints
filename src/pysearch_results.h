#include <Python.h>
#include "chemfp.h"

typedef struct {
    PyObject_HEAD
    int num_results;
    chemfp_search_result *results;
    PyObject *ids;
} SearchResults;

extern PyTypeObject chemfp_py_SearchResultsType;
