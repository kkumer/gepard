#define PY_SSIZE_T_CLEAN
#include <Python.h>

int main(int argc, char *argv[])
{
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    int i;
	double cff;
	double res[11];

    if (argc < 3) {
        fprintf(stderr,"Usage: kmcffs xB t Q2\n");
        return 1;
    }

    Py_Initialize();
    pName = PyUnicode_DecodeFSDefault("kmcffs");
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "getcffs");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 1);
            for (i = 0; i < argc - 1; ++i) {
				res[i] = atof(argv[i + 1]);
                pValue = PyFloat_FromDouble(atof(argv[i + 1]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
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
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "multiply");
        return 1;
    }
    if (Py_FinalizeEx() < 0) {
        return 120;
    }
    return 0;
}
