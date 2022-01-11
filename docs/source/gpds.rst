#####################
GPDs and form factors
#####################


Generic Gepard model
--------------------

GPDs and various form factors (CFFs for DVCS, TFFs for DVMP),
are all subclasses of the generic ``ParameterModel`` class. This means
that they all have a Python dictionary ``parameters`` with names of parameters and
their numerical values. For example, GPD model ``PWNormGPD`` (see the next
section) has among its parameters `ms2` which is squared sea quark
mass parameter controlling dipole :math:`t`-dependence:

.. code-block:: python

   >>> import gepard as g
   >>> gpd = g.PWNormGPD()
   >>> gpd.parameters['ms2']
   1.1

For more on manipulations with model parameters, see the 
:ref:`section on fitting<sec-fitting>`.

Apart from these "fitting" parameters, models usually have
additional attributes controlling their other properties, like
order and scheme of perturbative QCD, defining input scale :math:`Q_{0}^2`,
etc. These other parameters can be set and accessed via keyword arguments
or as attributes:

.. code-block:: python

   >>> gpd.p    # p=0 for LO, 1 for NLO
   0
   >>> gpd.scheme  # default factorization scheme
   'msbar'
   >>> gpd2 = g.PWNormGPD(p=1, q02=4)

Here we constructed second GPD model, `gpd`, which will be used
at NLO order, and considered to be defined at input scale
:math:`Q_{0}^2 = 4\, {\rm GeV}^2`,


GPDs
----

Generalized parton distributions (GPDs) are modelled in the
conformal moment :math:`j` space, with GPD of flavor :math:`a`,
:math:`H^{a}_{j}(\eta, t)`,
represented as ``H(eta, t)`` function of skewedness :math:`\eta` and
momentum transfer squared :math:`t`, which returns the two-dimensional
array of numbers with the shape :math:`npts \times 4`, where :math:`npts` is
the number of points on the Mellin-Barnes contour (which defines the
Compton form factors (CFFs)), and the four flavors are:

   - singlet quark
   - gluon
   - u_valence
   - d_valence

Valence here means "valence-like GPD" (see 
`hep-ph/0703179 <https://arXiv.org/abs/hep-ph/0703179>`_ for definition).

GPD :math:`E` is modeled in the same way, while GPDs :math:`\tilde{H}` and
:math:`\tilde{E}`, are not yet represented in Gepard on the GPD level
(corresponding CFFs are sometimes modeled directly, see below).

The conformal space GPDs :math:`H` and :math:`E`, are convoluted with
appropriate hard-scattering Wilson coefficients and evolution operator
to create CFFs (for DVCS) or transition form factors (TFFs, for DVMP).

The default conformal space GPD model in Gepard is ``PWNormGPD``,
which uses the SO(3) partial waves decomposition of GPDs, with
second and third partial wave proportional to the first one
(the only additional parameters are normalization of the waves).
This model is described in `arXiv:0904.0458 <https://arxiv.org/abs/0904.0458>`_.
The following table lists some more important attributes of this model.

.. list-table:: Attributes of PWNormGPD model
   :header-rows: 1

   * - Name
     - Description
     - Default value
   * - ``p``
     - pQCD order (0=LO, 1=NLO)
     - 0
   * - ``scheme``
     - pQCD factorization scheme ('msbar' or 'csbar')
     - 'msbar'
   * - ``nf``
     - number of active quark flavors
     - 4
   * - ``q02``
     - Initial :math:`Q_{0}^2` for GPD evolution
     - 4 (GeV^2)
   * - ``residual``
     - Residual :math:`t` dependence ('dipole' or 'exp')
     - 'dipole'


Here residual :math:`t` dependence is beyond the dependence
:math:`x^{-\alpha' t}` coming from the Regge trajectory.


CFFs
----

Compton form factors (CFFs) are required for calculation of DVCS observables.
There are presently four main classes of CFFs implemented in Gepard


.. _tab-CFF_classes:

.. list-table:: CFF classes
   :header-rows: 1

   * - Name of the class
     - Description
   * - ``MellinBarnesCFF``
     - Uses conformal-space GPDs as described in the section above
   * - ``DispersionFreePoleCFF``
     - Imaginary part of CFF is directly modelled, while real part
       is obtained by dispersion relations. :math:`\mathfrak{Re}\tilde{\mathcal{E}}` is given by freely
       parametrized pion pole.
   * - ``HybridCFF``
     - Combines ``MellinBarnesCFF`` for sea partons, with
       ``DispersionFreePoleCFF`` for valence quarks
   * - ``GoloskokovKrollCFF``
     - Model of Goloskokov and Kroll
       

Apart from ``GoloskokovKrollCFF`` which is completely fixed, all other CFFs
depend on parameters, either directly (``DispersionFreePoleCFF``) or via their
GPD model (``MellinBarnesCFF``). 

For CFFs which depend on GPD model, it is necessary to combine codes for
both CFF and GPD to get a working CFF model. This is done by creating a new
class like this:

.. code-block:: python

   >>> class MyCFF(g.PWNormGPD, g.MellinBarnesCFF):
   ...     pass
   >>> cff = MyCFF()

This is now a complete object and CFFs can be evaluated for some kinematics:

.. code-block:: python

   >>> cff.ReH(g.DataPoint({'xB': 0.1, 't': -0.3, 'Q2': 6}))
   13.44851


TFFs
----

Transition form factors (TFFs), are analogous to CFFs, but for DVMP,
and they in principle include, besides GPD, also a distribution
amplitude (DA) for the produced meson.
Presently, only the simplest TFF model is implemented, where DA is
given by its asymptotic form, while rest of the model is analogous

.. list-table:: TFF class
   :header-rows: 1

   * - Name of the class
     - Description
   * - ``MellinBarnesTFF``
     - Uses conformal-space GPDs and asymptotic DA


EFFs
----

To calculate DVCS observables, we also need elastic electromagnetic
form factors. There are three implementations in Gepard:

.. _tab-EFF_classes:

.. list-table:: EFF classes
   :header-rows: 1

   * - Name of the class
     - Description
   * - ``DipoleEFF``
     - dipole form of t-dependence
   * - ``KellyEFF``
     - EFFs as parametrized by J.J. Kelly, PRC 70 (2004) 068202
   * - ``ZeroEFF``
     - All EFFs are set to zero. Convenient for calculation of pure DVCS effects.


.. code-block:: python

   >>> eff = g.KellyEFF()
   >>> eff.F1(g.DataPoint({'t': 0}))  # Dirac form factor for proton
   1.0

