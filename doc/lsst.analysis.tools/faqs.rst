.. _analysis-tools-faqs:

FAQs
====

What Order are Configs Applied In?
----------------------------------

In analysis tools there are a number of ways to set config options which are applied in a certain order. If 
you get this order wrong you might be suprised by what your action is doing. In order from first applied to
last applied:

1. The action defaults
2. The config file in the obs package
3. The camera specific config in the obs package
4. The pipeline file configs
5. The command line --config and --config-file options

Due to this configs that are obs package specific should not be specified in the pipeline file because they 
overwrite the obs specific configs. For example the bands option, if that is set in the pipeline file it will 
overwrite the obs package which means that if you try to apply the same pipeline to HSC and DC2 data (for
example) it will not work for both as each survey has different band coverage.

How flag bands work
-------------------

The various flag selectors allow you to specify which bands the flags are applied in. If you want the selector
to be applied in the bands that the plot is being made in then the flag config option should be set to [].

Where is the line between python and YAML?
------------------------------------------

This can be a bit of a blurry line and in the end is up to the individual developer. A rough guideline is that 
if the code is reusable either for variants of the original action or useful for other actions then it belongs
in its own python class. If it is single use or a specific instance of a class then it should be done through
the YAML config. An example of this is the photometric repeatability metrics. Here there is a base action
called PhotometricRepeatability, defined in python, which is then configured for each individual application
in the pipeline YAML. 

.. code-block:: yaml

   atools.modelPhotRepGalSn5to10: PhotometricRepeatability
   atools.modelPhotRepGalSn5to10.fluxType: cModelFlux
   atools.modelPhotRepGalSn5to10.process.filterActions.perGroupStdevFiltered.selectors.extendedness.op: gt
   atools.modelPhotRepGalSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum: 5
   atools.modelPhotRepGalSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum: 10

The first line calls the action defined in python, after that the specific configuration is set in the
pipeline YAML. This is a good balance between defining the entire action in YAML and having duplicate actions 
defined in the python with only very mild differences.
