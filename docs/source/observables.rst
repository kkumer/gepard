#########################
Processes and observables
#########################


Processes
---------

Presently, only Deeply virtual Compton scattering (DVCS) process
can be considered fully implemented (using the so-called
BMK formulas). There is also an initial implementation of Deeply virtual
meson production (DVMP) of rho0 meson, valid for small Bjorken x
kinematics.

Datapoints should have an attribute ``process``, which is one of
the following:

.. _tab-processes:

.. table:: Processes
    :widths: auto

    +--------------------------------------+--------------------------+------------------+
    | Proces                               | DataPoint ``process``    | Generic process  |
    |                                      | value                    | class            |
    +======================================+==========================+==================+
    | :math:`e p \to e p \gamma`           | ``ep2epgamma``           | DVCS             |
    +--------------------------------------+--------------------------+------------------+
    | :math:`e n \to e n \gamma`           | ``en2engamma``           | DVCS             |
    +--------------------------------------+--------------------------+------------------+
    | :math:`\gamma^* p \to \gamma p`      | ``gammastarp2gammap``    | DVCS             |
    +--------------------------------------+--------------------------+------------------+
    | :math:`\gamma^* p \to \rho^{0} p`    | ``gammastarp2rho0p``     | DVMP             |
    +--------------------------------------+--------------------------+------------------+
    | :math:`e p \to e p X`                | ``dis``                  | DIS              |
    +--------------------------------------+--------------------------+------------------+

Formulas for generic theory class DVMP are implemented as
``g.DVMP``, while those for DVCS have several slightly different
versions:

.. _tab-BMK_formulas:

.. table:: Version of BMK formulas
    :widths: auto

    +--------------------+------------------------------------------------------------------+
    | Theory class       | description                                                      |
    +====================+==================================================================+
    | ``g.BMK``          | old original BMK formulas from 2005                              |
    +--------------------+------------------------------------------------------------------+
    | ``g.hotfixedBMK``  | improved formulas that are better for JLab kinematics            |
    +--------------------+------------------------------------------------------------------+
    | ``g.BM10ex``       | best formulas from 2010                                          |
    +--------------------+------------------------------------------------------------------+
    | ``g.BM10``         | like ``BM10ex``, but with some Q2-suppresed terms removed        |
    +--------------------+------------------------------------------------------------------+
    | ``g.BM10tw2``      | like ``BM10`` but with higher twists set to zero                 |
    +--------------------+------------------------------------------------------------------+


.. _tab-observables:

Observables
-----------

Following observables are implemented in Gepard. Identifiers
from the table can be used either
as methods of theory objects, or as ``observable`` attributes of datapoints.

.. table:: DVCS observables
    :widths: auto

    +----------------+------------------------------------------------------------------+
    | Identifier     | Description                                                      |
    +================+==================================================================+
    | ``XS``         | cross-section for leptoproduction of real photon                 |
    +----------------+------------------------------------------------------------------+
    | ``XGAMMA``     | cross-section for production of real photon by virtual one       |
    +----------------+------------------------------------------------------------------+
    | ``XUU``        | beam spin sum a.k.a helicity independent XS                      |
    +----------------+------------------------------------------------------------------+
    | ``XUUw``       | XUU weighted by BH propagator                                    |
    +----------------+------------------------------------------------------------------+
    | ``XLU``        | beam spin difference a.k.a helicity dependent XS                 |
    +----------------+------------------------------------------------------------------+
    | ``XLUw``       | XLU weighted by BH propagator                                    |
    +----------------+------------------------------------------------------------------+
    | ``XCLU``       | beam charge-spin difference (COMPASS)                            |
    +----------------+------------------------------------------------------------------+
    | ``XCUU``       | beam charge-spin sum (COMPASS)                                   |
    +----------------+------------------------------------------------------------------+
    | ``AC``         | beam charge asymmetry                                            |
    +----------------+------------------------------------------------------------------+
    | ``ALU``        | beam spin asymmetry                                              |
    +----------------+------------------------------------------------------------------+
    | ``ALUI``       | beam spin asymmetry, interference part                           |
    +----------------+------------------------------------------------------------------+
    | ``ALUDVCS``    | beam spin asymmetry, DVCS part                                   |
    +----------------+------------------------------------------------------------------+
    | ``AUL``        | longitudinal target spin asymmetry                               |
    +----------------+------------------------------------------------------------------+
    | ``AUT``        | transversal target spin asymmetry                                |
    +----------------+------------------------------------------------------------------+
    | ``AUTI``       | transversal target spin asymmetry, interference part             |
    +----------------+------------------------------------------------------------------+
    | ``AUTDVCS``    | transversal target spin asymmetry, DVCS part                     |
    +----------------+------------------------------------------------------------------+
    | ``BTSA``       | beam (longitudinal) target double spin asymmetry                 |
    +----------------+------------------------------------------------------------------+
    | ``ALTI``       | beam transversal target double spin asymmetry, interference part |
    +----------------+------------------------------------------------------------------+
    | ``ALTBHDVCS``  | beam transversal target double spin asymmetry, BH-DVCS part      |
    +----------------+------------------------------------------------------------------+

Many of these observables can be evaluated both as differential in azimuthal
angle :math:`\phi` (if the ``DataPoint`` argument has an attribute ``phi``),
or as "harmonic", i. e., as Fourier integral over :math:`\phi` (if the
``DataPoint`` argument has attribute ``FTn``).
Similary, ``XGAMMA`` will be evaluated as differential in :math:`t` if
``DataPoint`` has attribute ``t``, and as integrated over :math:`t` if
it doesn't.

.. table:: DVMP observables
    :widths: auto

    +--------------------+------------------------------------------------------------------+
    | Name               | Description                                                      |
    +====================+==================================================================+
    | ``XGAMMA``         | cross-section for production of meson by longit. virtual photon  |
    +--------------------+------------------------------------------------------------------+

Choice whether DVCS or DVMP ``XGAMMA`` will be evaluated is dependent
on the value of ``pt.process``.
