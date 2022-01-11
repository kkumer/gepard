###################
Building the theory
###################

``Theory`` object is in some sense a collection of formulas
which implement that what you want to calculate.

Several complete theories (for DVCS only, corresponding to the
so-called "KM" models) are already available in Gepard
and can be imported and used like this

.. code-block:: python

   >>> from gepard.fits import th_KM15

but when user wants to go beyond that, they should build the
``Theory`` object themselves.

For convenience, formulas are separated into building blocks, so
that user can make choice among different options. Procedure
for building a new theory goes roughly as follows.

First you choose the GPD model. Presently, the only choices are ``g.PWNormGPD`` and ``g.TestGPD``.

Then think whether you are interested in DVCS, DVMP, DIS, or some combination of these
processes? For DVCS you then need to choose the
CFF model. Available choices are listed in table :ref:`CFF classes<tab-CFF_classes>`. 
You also need
to chose model for elastic form factors, among those listed in table :ref:`EFF classes<tab-EFF_classes>`.
For DVMP, the only available choice for TFFs is ``g.MellinBarnesTFF``.
Finally, you chose expressions for observables. For DVCS, choice is among different
versions of BMK formulas, listed in table :ref:`BMK formulas<tab-BMK_formulas>`, while for DVMP,
there is only ``g.DVMP``. For DIS F2, there is ``g.DIS``.


Then, from all these building blocks, you build a complete theory like this:


.. code-block:: python

   >>> import gepard as g
   >>> class MyTheory(g.PWNormGPD, g.MellinBarnesCFF, g.KellyEFF, g.BM10tw2, g.MellinBarnesTFF, g.DVMP, g.DIS):
   ...     pass
   >>> th = MyTheory()


This is an example of the "maksimal" theory, which includes everything that Gepard can calculate.
There is almost no performance hit by including what you don't need. Also, the order of
classes in the parentheses should not be important, because care has been taken that there
are no name clashes between things implemented by different classes.
