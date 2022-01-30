########
Tutorial
########


After importing Gepard package, if one of the theories available
in ``gepard.fits`` are not enough, and we want to make a new
GPD fit, we can construct a new ``theory`` class 
by combining code for

1. GPDs (we take SO(3) partial-waves conformal space model ``PWNormGPD``), 
2. for CFFs ( we take ``MellinBarnesCFF`` which combines conformal space GPDs with appropriate hard-scattering coefficients), and
3. for DVCS observables (``BMK``) 
4. for DIS F2 (good idea, so that forward limit of GPDs describes DIS data)

.. plot::
   :context:
   :include-source:

   >>> import gepard as g
   >>> class MyTheory(g.PWNormGPD, g.MellinBarnesCFF, g.DIS, g.BMK):
   ...     pass


We then construct an instance of this new class, where we
decide that we want to work at NLO (``p = 1``):

.. plot::
   :context:
   :include-source:

   >>> th = MyTheory(p=1)

.. note::

   It is somewhat of a convention, which we follow in this documentation,
   that names of ``Theory`` instances always start with ``th``, while
   the names of ``DataPoint`` instances always start with ``pt``.


We would now like to confront this theory with experimental data on DIS
and DVCS. For demonstration purposes, we take just two H1 datasets:

.. plot::
   :context:
   :include-source:

   >>> DISpoints = g.dset[206]
   >>> DVCSpoints = g.dset[39]

To see in more detail what are these datasets, we can use utility
function ``describe_data``:

.. plot::
   :context:
   :include-source:

   >>> g.describe_data(DISpoints+DVCSpoints)
   npt x obs     collab  FTn    id  ref.        
   ----------------------------------------------
   10 x DISF2   H1      N/A    206 Nucl.Phys.B470(96)3
    8 x XGAMMA  H1      N/A    39  hep-ex/0505061
   ----------------------------------------------
   TOTAL = 18

Here we see that dataset ``id=206`` contains 10
measurements of DIS :math:`F_2`, while ``id=39`` contains
8 measurements of pure DVCS cross-section, differential in :math:`t`.
(``FTn`` column is not relevant here.). We also see the
reference to original literature where measurements were published.

To see if our theory describes this data we can either calculate
total :math:`\chi^2`:

.. plot::
   :context:
   :include-source:

   >>> th.chisq(DVCSpoints)
   85.21

or we can plot theory line against data points like this:

.. plot::
   :context:
   :include-source:

   >>> import gepard.plots
   >>> gepard.plots.jbod(points=DVCSpoints, lines=th).show()


This is obviously bad, so let us fit the parameters of the
theory to this data. For this, we construct the ``MinuitFitter`` object,
release some of the model parameters (overal normalization ``ns``,
residual :math:`t`-dependence parameter ``ms2``, and normalization
of the second partial wave ``secs``, all for sea quarks):

.. plot::
   :context: close-figs
   :include-source:

   >>> f = g.MinuitFitter(DISpoints+DVCSpoints, th)
   >>> f.release_parameters('ns', 'ms2', 'secs')
   >>> f.fit()

After fitting is done, we print the resulting values and uncertainties of fitting parameters:

.. plot::
   :context: close-figs
   :include-source:

   >>> th.print_parameters()
   ns    =    0.17 +- 0.01
   ms2   =    0.93 +- 0.10
   secs  =    0.18 +- 0.03

Theory now describes the data fine, as one can see from :math:`\chi^2`
value:

.. plot::
   :context: close-figs
   :include-source:

   >>> th.chisq(DISpoints+DVCSpoints)
   6.33

and, visually, from the plot:


.. plot::
   :context: close-figs
   :include-source:

   >>> gepard.plots.jbod(points=DVCSpoints, lines=th).show()


Finally, one could calculate and then plot some particular CFF, like this:


.. plot::
   :context: close-figs
   :include-source:

   >>> import numpy as np
   >>> import matplotlib.pyplot as plt
   >>> xis = np.linspace(0.001, 0.1)
   >>> ims = []
   >>> res = []
   >>> for xi in xis:
   ...     pt = g.DataPoint(xi=xi, t=-0.2, Q2=4)
   ...     ims.append(xi*th.ImH(pt))
   ...     res.append(xi*th.ReH(pt))
   >>> plt.plot(xis, ims, label='Im(H)')  # doctest: +SKIP
   >>> plt.plot(xis, res, label='Re(H)')  # doctest: +SKIP
   >>> plt.xlabel(r'$\xi$', fontsize=14)  # doctest: +SKIP
   >>> plt.ylabel('ImH', fontsize=14)  # doctest: +SKIP
   >>> plt.legend()  # doctest: +SKIP


