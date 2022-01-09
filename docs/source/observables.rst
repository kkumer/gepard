#########################
Processes and observables
#########################

Presently, only Deeply virtual Compton scattering (DVCS) process
can be considered fully implemented (using the so-called
BMK formulas). There is also an implementation of Deeply virtual
meson production (DVMP) of rho meson, valid for small Bjorken x
kinematics.


.. _tab-observables:

Observables
-----------

Following observables are implemented in gepard. They can be used
as methods of theory objects, or ``yaxis`` attributes of datapoints.

.. list-table:: Observables
   :header-rows: 1

   * - Name
     - Description
   * - XS
     - Cross-section
   * - BSS
     - beam spin sum a.k.a helicity independent XS
   * - BSSw
     - BSS weighted by BH propagator
   * - BSD
     - beam spin difference a.k.a helicity dependent XS
   * - BSDw
     - BSD weighted by BH propagator
   * - BCA
     - beam charge asymmetry
   * - BSA
     - beam spin asymmetry
   * - ALUI
     - beam spin asymmetry, interference part
   * - ALUDVCS
     - beam spin asymmetry, DVCS part
   * - TSA
     - (longitudinal) target spin asymmetry
   * - AUTI
     - transversal target spin asymmetry, interference part
   * - AUTDVCS
     - transversal target spin asymmetry, DVCS part
   * - BTSA
     - beam (longitudinal) target double spin asymmetry
   * - ALTI
     - beam transversal target double spin asymmetry, interference part
   * - ALTBHDVCS
     - beam transversal target double spin asymmetry, BH-DVCS part


