#########################
Processes and observables
#########################


Processes
---------

Presently, only Deeply virtual Compton scattering (DVCS) process
can be considered fully implemented (using the so-called
BMK formulas). There is also an initial implementation of Deeply virtual
meson production (DVMP) of rho meson, valid for small Bjorken x
kinematics.

Datapoints should have an attribute ``process``, which is one of
the following:

.. _tab-processes:

.. list-table:: Processes
   :header-rows: 1

   * - Proces
     - DataPoint ``process`` value
     - Generic theory class
   * - :math:`e p \to e p \gamma`
     - ``ep2epgamma``
     - DVCS
   * - :math:`e n \to e n \gamma`
     - ``en2engamma``
     - DVCS
   * - :math:`\gamma^* p \to \gamma p`
     - ``gammastarp2gammap``
     - DVCS
   * - :math:`\gamma^* p \to \rho^{0} p`
     - ``gammastarp2rho0p``
     - DVMP


Formulas for generic theory class DVMP are implemented as
``g.DVMP``, while those for DVCS have several slightly different
versions:

.. _tab-BMK_formulas:

.. list-table:: Version of BMK formulas
   :header-rows: 1

   * - Theory class
     - description
   * - ``g.BMK``
     - old original BMK formulas from 2005
   * - ``g.hotfixedBMK``
     - improved formulas that are better for JLab kinematics
   * - ``g.BM10ex``
     - best formulas from 2010
   * - ``g.BM10``
     - like ``BM10ex``, but with some Q2-suppresed terms removed
   * - ``g.BM10tw2``
     - like ``BM10`` but with higher twists set to zero




.. _tab-observables:

Observables
-----------

Following observables are implemented in gepard. They can be used
as methods of theory objects, or ``yaxis`` attributes of datapoints.

.. list-table:: DVCS observables
   :header-rows: 1

   * - Identifier
     - Description
   * - ``XS``
     - cross-section for leptoproduction of real photon
   * - ``XDVCS``
     - cross-section for production of real photon by virtual one
   * - ``BSS``
     - beam spin sum a.k.a helicity independent XS
   * - ``BSSw``
     - BSS weighted by BH propagator
   * - ``BSD``
     - beam spin difference a.k.a helicity dependent XS
   * - ``BSDw``
     - BSD weighted by BH propagator
   * - ``BCA``
     - beam charge asymmetry
   * - ``BSA``
     - beam spin asymmetry
   * - ``ALUI``
     - beam spin asymmetry, interference part
   * - ``ALUDVCS``
     - beam spin asymmetry, DVCS part
   * - ``TSA``
     - (longitudinal) target spin asymmetry
   * - ``AUTI``
     - transversal target spin asymmetry, interference part
   * - ``AUTDVCS``
     - transversal target spin asymmetry, DVCS part
   * - ``BTSA``
     - beam (longitudinal) target double spin asymmetry
   * - ``ALTI``
     - beam transversal target double spin asymmetry, interference part
   * - ``ALTBHDVCS``
     - beam transversal target double spin asymmetry, BH-DVCS part


.. list-table:: DVMP observables
   :header-rows: 1

   * - Name
     - Description
   * - ``X``
     - cross-section for production of meson by virtual photon
