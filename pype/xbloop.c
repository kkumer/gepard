#define PY_SSIZE_T_CLEAN
#include <Python.h>

int main(int argc, char *argv[])
{
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    int i, k, nxB;
	double xB, xBmin, xBmax, step, t, Q2, cff;
	double res[11];

    t = -0.3;
    Q2 = 4.0;
    xBmin = 0.01;
    xBmax = 0.9;
    nxB = 40;
    step = (xBmax-xBmin)/(nxB-1.);

    if (argc != 1) {
        fprintf(stderr,"Usage: xbloop\n");
        return 1;
    }

    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    pName = PyUnicode_DecodeFSDefault("kmcffs");
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    /* header */
      printf("%s", "   xB        t         Q2      ImH      ReH      ImE      ReE      ImHt     ReHt    ImEt      ReEt\n");

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "getcffs");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            for (k = 0; k < nxB; ++k) {
            /* arguments for Python function */
                pArgs = PyTuple_New(3);
                xB = xBmin + k*step;
                pValue = PyFloat_FromDouble(xB);
                PyTuple_SetItem(pArgs, 0, pValue);
                pValue = PyFloat_FromDouble(t);
                PyTuple_SetItem(pArgs, 1, pValue);
                pValue = PyFloat_FromDouble(Q2);
                PyTuple_SetItem(pArgs, 2, pValue);
                /* calling Python function */
                pValue = PyObject_CallObject(pFunc, pArgs);
                Py_DECREF(pArgs);
                if (pValue != NULL) {
                    PyObject *iter = PyObject_GetIter(pValue);
                    if (!iter) {
                        fprintf(stderr, "Must return list/tuple");
                    }
                    i = 0;
                    while (1) {
                        PyObject *next = PyIter_Next(iter);
                          if (!next) {
                            /* nothing left in the iterator  */
                            break;
                          }

                          if (!PyFloat_Check(next)) {
                            fprintf(stderr,"Non-float returned!\n");
                            break;
                          }
                          cff = PyFloat_AsDouble(next);
                          res[i] = cff;
                          i++;
                    }
                    for (i = 0; i < 11; i++){
                        printf("%f ", res[i]);
                    }
                    printf("%s", "\n");
                    Py_DECREF(pValue);
                }
                else {
                    Py_DECREF(pFunc);
                    Py_DECREF(pModule);
                    PyErr_Print();
                    fprintf(stderr,"Call failed\n");
                    return 1;
                }
          }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", "getcffs");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load Python module \"%s\"\n", "kmcffs");
        return 1;
    }
    if (Py_FinalizeEx() < 0) {
        return 120;
    }
    return 0;
}
