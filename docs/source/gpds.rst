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
   >>> th = g.PWNormGPD()
   >>> th.parameters['ms2']
   1.1

For more on manipulations with model parameters, see the 
:ref:`section on fitting<sec-fitting>`.

Apart from these "fitting" parameters, models usually have
additional attributes controlling their other properties, like
order and scheme of perturbative QCD, defining input scale :math:`Q_{0}^2`,
etc. These other parameters can be set and accessed via keyword arguments
or as attributes:

.. code-block:: python

   >>> th.p    # p=0 for LO, 1 for NLO
   0
   >>> th.scheme  # default factorization scheme
   'msbar'
   >>> th2 = g.PWNormGPD(p=1, Q02=2)

Here we constructed second GPD model, ``th2``, which will be used
at NLO order, and is considered to be defined at lower input scale
:math:`Q_{0}^2 = 2\, {\rm GeV}^2`.
(Most of the Gepard models have :math:`Q_{0}^2 = 4\, {\rm GeV}^2`
as a default input scale.)


GPDs
----

Generalized parton distributions (GPDs) are modelled in the
conformal moment :math:`j` space, with GPD of flavor :math:`a`,
:math:`H^{a}_{j}(\eta, t)`,
implemented as ``H(eta, t)`` Python function of skewedness :math:`\eta` and
momentum transfer squared :math:`t`, which returns the two-dimensional
array of numbers with the shape :math:`npts \times 4`, where :math:`npts` is
the number of points on the Mellin-Barnes contour (the complex-space
contour which defines the Compton form factors (CFFs)), and the 4 flavors are:

   - singlet quark
   - gluon
   - u_valence
   - d_valence

Valence here means "valence-like GPD" (see 
`hep-ph/0703179 <https://arXiv.org/abs/hep-ph/0703179>`_ for definition).

GPD :math:`E` is modeled in the same way, while GPDs :math:`\tilde{H}` and
:math:`\tilde{E}`, are not yet represented in Gepard on the GPD level.
(In some Gepard theories corresponding CFFs :math:`\tilde{\tilde{H}}`,
and :math:`\tilde{\tilde{E}}` are modeled directly.)

The conformal space GPDs :math:`H` and :math:`E`, are convoluted with
appropriate hard-scattering Wilson coefficients and evolution operator
to create CFFs (for DVCS) or transition form factors (TFFs, for DVMP).

The default conformal space GPD model in Gepard is ``PWNormGPD``,
which uses the SO(3) partial waves decomposition of GPDs, with
second and third partial wave proportional to the first one
(the only additional parameters of subleading partial waves are their normalizations).
This model is described in `arXiv:0904.0458 <https://arxiv.org/abs/0904.0458>`_.
The following table lists some more important attributes of this model.

.. table:: Attributes of PWNormGPD model
    :widths: auto

    +------------------+--------------------------------------------------------+------------------+
    | Name             | Description                                            | Default value    |
    +==================+========================================================+==================+
    | ``p``            | pQCD order (0=LO, 1=NLO)                               | 0                |
    +------------------+--------------------------------------------------------+------------------+
    | ``scheme``       | pQCD factorization scheme ('msbar' or 'csbar')         | 'msbar'          |
    +------------------+--------------------------------------------------------+------------------+
    | ``nf``           | number of active quark flavors                         | 4                |
    +------------------+--------------------------------------------------------+------------------+
    | ``Q02``          | Initial :math:`Q_{0}^2` for GPD evolution              | 4 (GeV^2)        |
    +------------------+--------------------------------------------------------+------------------+
    | ``residual``     | Residual :math:`t` dependence ('dipole' or 'exp')      | 'dipole'         |
    +------------------+--------------------------------------------------------+------------------+

Here the residual :math:`t` dependence is additional to the dependence
:math:`x^{-\alpha' t}` coming from the Regge trajectory.
Parameters of this model are described in the above paper.

For a given j-space model, you can evaluate also standard GPDs
in the x-space (presently only for :math:`\eta=0`
or :math:`\eta=x`), using method ``th.Hx``,
which returns the triplet [singlet/sea quark, gluon,
non-singlet/valence quark] GPDs (the third one is not implemented yet
and is set to zero):

.. code-block:: python

   >>> pt = g.DataPoint(x=0.01, eta=0, t=0, Q2=8)
   >>> th.Hx(pt)   # should be equal to PDFs
   array([258.04908329,  25.19998599,   0.        ])
   >>> pt = g.DataPoint(x=0.01, eta=0.01, t=-0.2, Q2=4)
   >>> th.Hx(pt)
   array([267.91985614,   2.13496561,   0.        ])
   >>> pt.eta = 0.3
   >>> th.Hx(pt)   # doesn't work yet for arbitrary eta!
   Traceback (most recent call last):
   ...
   Exception: eta has to be either 0 or equal to x


CFFs
----

Compton form factors (CFFs) are required for calculation of DVCS observables.
There are presently four main classes of CFFs implemented in Gepard


.. _tab-CFF_classes:

.. table:: CFF classes
    :widths: auto

    +--------------------------+------------------------------------------------------------------+
    | Name of the class        | Description                                                      |
    +==========================+==================================================================+
    | ``MellinBarnesCFF``      | Uses conformal-space GPDs as described in the section above      |
    +--------------------------+------------------------------------------------------------------+
    | ``DispersionFreePoleCFF``| Imaginary part of CFF is directly modelled, while real part is   |
    |                          | obtained by dispersion relations.                                |
    |                          | :math:`\mathfrak{Re}\tilde{\mathcal{E}}` is given by freely      |
    |                          | parametrized pion pole.                                          |
    +--------------------------+------------------------------------------------------------------+
    | ``HybridCFF``            | Combines ``MellinBarnesCFF`` for sea partons, with               |
    |                          | ``DispersionFreePoleCFF`` for valence quarks                     |
    +--------------------------+------------------------------------------------------------------+
    | ``GoloskokovKrollCFF``   | ``GoloskokovKrollCFF`` Model of Goloskokov and Kroll             |
    +--------------------------+------------------------------------------------------------------+
       

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

   >>> cff.ReH(g.DataPoint(xB=0.1, t=-0.3, Q2=6))
   13.44851


TFFs
----

Transition form factors (TFFs), are analogous to CFFs, but for DVMP,
and they in principle include, besides GPD, also a distribution
amplitude (DA) for the produced meson.
Presently, only the simplest TFF model is implemented, where DA is
given by its asymptotic form, while rest of the model is analogous

.. table:: TFF class
    :widths: auto

    +--------------------------+--------------------------------------------------------------+
    | Name of the class        | Description                                                  |
    +==========================+==============================================================+
    | ``MellinBarnesTFF``      | Uses conformal-space GPDs and asymptotic DA                  |
    +--------------------------+--------------------------------------------------------------+


EFFs
----

To calculate DVCS observables, we also need elastic electromagnetic
form factors. There are three implementations in Gepard:

.. _tab-EFF_classes:

.. table:: EFF classes
    :widths: auto

    +--------------------+--------------------------------------------------------------+
    | Name of the class  | Description                                                  |
    +====================+==============================================================+
    | ``DipoleEFF``      | dipole form of t-dependence                                  |
    +--------------------+--------------------------------------------------------------+
    | ``KellyEFF``       | EFFs as parametrized by J.J. Kelly, PRC 70 (2004) 068202     |
    +--------------------+--------------------------------------------------------------+
    | ``ZeroEFF``        | All EFFs are set to zero. Convenient for calculation of pure |
    |                    | DVCS effects.                                                |
    +--------------------+--------------------------------------------------------------+


.. code-block:: python

   >>> eff = g.KellyEFF()
   >>> eff.F1(g.DataPoint(t=0))  # Dirac form factor for proton
   1.0

