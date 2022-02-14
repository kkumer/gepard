.. _sec-fitting:

######################
Fitting theory to data
######################


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
   8.411
   >>> f.print_parameters()
   ns    =    1.46 +- 0.45
   ms2   =    0.93 +- 0.07
   secs  =   -0.32 +- 0.01


After successful fit of theory object ``th``, user can access parameter uncertainties as
``th.parameters_errors`` dictionary, and full covariance matrix (inverse of the
chi-square Hessian matrix) as ``th.covariance`` dictionary.

Covariance matrix can then be used to propagate uncertainty to prediction
of observables, like this:


   >>> th.predict(pts[0], uncertainty=True)
   (13.25, 1.62)


where parameter-dependent form factors (such as CFFs) can also be "predicted",
i. e.,  calculated together with their uncertainty:

  >>> pt = g.DataPoint(xB=0.01, t=-0.2, Q2=10)
  >>> th.predict(pt, observable='ImH', uncertainty=True)
  (273.2, 25.8)


.. note::

   These uncertainties are "naive". They ignore unknown part of uncertainty due
   to the rigidity of the model. It is likely that true uncertainty is
   always significantly larger.


Fitting using ``f.fit()`` is simple but limited. For better control over fitting procedure
and determination of parameter uncertainties, user should use many functionalities
of the ``iminuit`` package, and call directly its functions like ``f.minuit.migrad``,
``f.minuit.minos`` etc.
For details one should consult  `documentation <https://iminuit.readthedocs.io/en/stable/>`_
of this package.
It is important to keep in mind that although parameter values are automatically kept in sync between
``iminuit`` and ``Gepard``, parameter uncertainties and covariances are *not*. So, in 
order to correctly use uncertainties in ``Gepard``, one should
first syncronize covariances using special function ``f.covsync``. For example, simple
fitting with ``f.fit()`` is equivalent to first doing ``f.minuit.migrad()`` and then
``f.covsync()``.

