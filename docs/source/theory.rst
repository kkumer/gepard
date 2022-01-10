###################
Building the theory
###################

``Theory`` object is in some sense a collection of formulas
which implement everything that you want to calculate.

For convenience, formulas are separated into building blocks that
can in principle be used separately

If you want to make a new model, beyond what is already available in
Gepard you can work like this:


.. code-block:: python

   >>> import gepard as g
   >>> class MyTheory(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):
   ...     pass
   ...
   ...     def ImE(pt):
   ...         return 0   
   >>> th = MyTheory()


