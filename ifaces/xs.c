#include <Python.h>
#include "ifaces.h"

/* FIXME: pValue not DECREF'd ! */

double xs(int id, int Q, int lam, double Ee, double Ep, double xB, double Q2, double t, double phi)
{
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue, *pResult;
    int i;
    double kin[9];
    double res;

    kin[0] = id; 
    kin[1] = Q;
    kin[2] = lam;
    kin[3] = Ee;
    kin[4] = Ep;
    kin[5] =  xB;
    kin[6] =  Q2;
    kin[7] =  t;
    kin[8] =  phi;
    res = 0.;

    Py_Initialize();
    pName = PyString_FromString("xs");
    /* name of the module */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "XSunp");
        /* callable function from the module */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(9);
            for (i = 0; i < 9; ++i) {
                pValue = PyFloat_FromDouble(kin[i]);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            /* Calling Python function now */
            pResult = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pResult != NULL) {
                res = PyFloat_AsDouble(pResult);
                /* printf("Result of call: %8.3e\n", res); */
                Py_DECREF(pResult);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", "XSunp");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "grid");
        return 1;
    }
    /* Py_Finalize(); */
    /* Py_Finalization is not allowed due to 
     * this bug: 
     * https://bugs.launchpad.net/ubuntu/+source/python-numpy/+bug/184920
     */
    return res;
}
