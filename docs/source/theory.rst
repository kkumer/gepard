###################
Building the theory
###################

``Theory`` object consists of collection of formulas
which implement everything that you want to calculate.

...

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


