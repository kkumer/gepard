#######
Fitting
#######


At the moment, only the standard least-squares fitting is implemented,
and it uses the `iminuit <https://iminuit.readthedocs.io/en/stable/>`_
Python interface for the Minuit2 C++ library.

Gepard user creates a ``MinuitFitter`` object, where first argument
is a dataset (collection of ``DataPoint`` objects), and second
argument is a ``Theory`` object. Then ``fit()`` method of this
``MinuitFitter`` objects starts the fit (by calling ``migrad``
optimizer of Minuit). Before this, user should release some
parameters of the ``Theory`` because all parameters are 
by default fixed at the moment of the creation of ``Theory``
object.

So, the minimal example of fitting is:


.. code-block:: python

   >>> import gepard as g
   >>> pts = g.dset[36]     # 12 H1 DVCS measurements
   >>> class FitTest(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):
   ...     pass
   >>> th = FitTest()
   >>> f = g.MinuitFitter(pts, th)
   >>> f.release_parameters('ns', 'ms2', 'secs')
   >>> f.fit()


Final values of chi-square and of parameters are available as

   >>> f.minuit.fval
   8.411290505223722
   >>> f.print_parameters()
   ns    =    1.455 +- 0.448
   ms2   =    0.927 +- 0.067
   secs  =   -0.315 +- 0.005


Using ``f.minuit`` user can directly access all the functionalities of the ``iminuit``,
and should consult its `documentation <https://iminuit.readthedocs.io/en/stable/>`_

After successful fit of theory object ``th``, user can access parameter uncertainties as
``th.parameters_errors`` dictionary, and full covariance matrix (inverse of the
chi-square Hessian matrix) as ``th.covariance`` dictionary.

Covariance matrix can then be used to propagate uncertainty to prediction
of observables, like this:


   >>> th.predict(pts[0], uncertainty=True)
   (13.253860520202533, 1.616903222077874)


where parameter-dependent form factors (such as CFFs) can also be "predicted",
i. e.,  calculated together with their uncertainty:

  >>> pt = g.DataPoint({'xB': 0.01, 't': -0.2, 'Q2':10})
  >>> th.predict(pt, observable='ImH', uncertainty=True)
  (273.1888807090146, 25.75226315614056)


.. note::

   These uncertainties are "naive". They ignore unknown part of uncertainty which
   results from the rigidity of the model. It is likely that true uncertainty is
   always significantly larger.



