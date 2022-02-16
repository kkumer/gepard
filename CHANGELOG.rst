Changelog
#########

v0.9.10 (2022-02-16)
--------------------

* Fix pip dependencies

v0.9.9 (2022-02-15 )
-------------------

* First public Python package release
* NLO DVEM works now (rho0 only)
* Removed (temporarily) NNLO and neural net stuff


v0.9.8 (2021-12-06)
-------------------

* First beta with proper distributable Python package structure
* Added GNU AGPL-3.0 license.


v0.9.7 (2015-10-07)
-------------------

* ADD: Merged pype into gepard
* CHG: Moved to git and gitlab.com


v0.9.6 (2013-08-07)
-------------------

* ADD: DVEM coefficient functions (just starting, DVEM not usable yet)
* ADD: GPD E
* FIX: All test for "big" NPB paper (2008) fixed


v0.9.5 (2010-06-25)
-------------------

* ADD: Evaluation of gamma function is also controlled by
	       SPEED, saving significant time in fitting
* ADD: Compiles also with Intel ifort, NAG f95 and GNU gfortran compilers,
	       including OpenMP parallelization.
* ADD: python extension pygepard.so created
* ADD: third SO(3) partial wave to 'FIT' model
* ADD: Evolution of GPDs at eta=0 and eta=x
* ADD: Possibility of eta-dependent ansatz
* ADD: Mathematica interface completely redone so both
	       calculation of CFFs and fitting is possible from
		   within Mathematica
* ADD: Evolved Wilson coefficients can be now calculated in
	       the initialization phase which leads to much faster
		   execution
* FIX: comparison to PEGASUS is OK at NNLO as well
* FIX: SPEED=3 should be used instead of SPEED=4 now


v0.9.3 (2006-12-17)
-------------------

* ADD: non-singlet evolution
* ADD: non-singlet NNLO
* ADD: MSbar non-diagonal evolution
* ADD: intelligent choice of y-range for drawing of PDF(x)
* ADD: many example programs in ex subdirectory which produce
	       data for plots in a big paper
* ADD: automatic plotting and psfragging of plots
* FIX: dependence on factorization scale corrected now
* BUG: SPEED=4 is much more unreliable than before
* FIX: charge-factor bugs (it was not NF-dependent, and was
	       altogher missing from DVCS)
* ADD: PDF(x), resulting from fit, are also drawn now
* ADD: Additional experimental data points for fitting


v0.9.2  (2006-08-03)
--------------------

* Fixed F2(P) bug that forced P=1 in DIS case


v0.9.1  (2006-07-26)
--------------------

* Added possibility of Gaussian integration on 2**ACC points
	  with new parameter ACC=1..6.
* Added test routines 'houches' and 'accuracy' for checking
	  accuracy w.r.t. Les Houches benchmark and w.r.t
	  referent most precise gepard result


v0.9  (2006-07-23)
------------------

* Wilson coefficients calculated in initialization phase and
	  together with anomalous dimensions put into common block
* Integration moved from slatec to 8-point Gaussian
* Adacf part integrated into source
* PGPLOT support added
* Compiles cleanly on both Linux and Cygwin
* Additional user-friendly documentation in gepard.pdf


v0.88  (2006-07-02)
-------------------

* Initial version
* Imported into subversion versioning system
* Target gepard.exe is broken (so calling from Mathematica
		doesn't work) but "pure Fortran" targets work and produce
		data equivalent with Figs 1 and 2 of Letter


v0.82 (2006-04-06)
------------------

* Numerical integration improved by subdivision of path:
* LAM..-1/2..1/2..LAM


v0.8 (2006-03-24)
-----------------

* Some old release


v0.7 (2006-03-20)
-----------------

* Some old release
