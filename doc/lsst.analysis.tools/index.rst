.. py:currentmodule:: lsst.analysis.tools

.. _lsst.analysis.tools:

###################
lsst.analysis.tools
###################

.. Paragraph that describes what this Python module does and links to related modules and frameworks.

.. .. _lsst.analysis.tools-using:

Using lsst.analysis.tools
=========================

.. toctree linking to topics related to using the module's APIs.

.. .. toctree::
..    :maxdepth: 1

Getting Started Guide
=====================

What is analysis_tools?
-----------------------

Analysis tools is a package for the creation of plots and metrics. It allows
them to be created from the same code to ensure that they are consistent
and repeatable.

It has a lot of very powerful and flexible functionality but can be a little
bit overwhelming at first. We are going to cover getting the stack set up 
and cloning analysis_tools. Then an overview of the package and follow that 
with walking through a series of examples of increasing complexity.

Analysis tools is designed to work with any sort of keyed data but to make it 
more intuitive initially we’ll talk about tables and column names.

Setting Up the Package and Getting Started With The Stack
---------------------------------------------------------

To set up the stack you need to source your version of:

``source /opt/lsst/software/stack/loadLSST.bash``

and then setup lsst_distrib.

``setup lsst_distrib``

If you don't have a local version of analysis_tools but you want to contribute or change things
then you can git clone the package from https://github.com/lsst/analysis_tools.

``git clone git@github.com:lsst/analysis_tools.git``

Once the package has been closed it needs to be setup and scons needs to be run.

``setup -j -r repos/analyis_tools``
``scons repos/analysis_tools``

More details can be found here:
https://pipelines.lsst.io/install/package-development.html?highlight=github#un-set-up-the-development-package
You can check what is set up using this command:

``eups list -s``

hopefully this will show your local version of analysis_tools.

Package Layout
==============

There are a bunch of files in analysis_tools but we are going to focus on two directories, 
``python/lsst/analysis/tools/`` and ``pipelines``, which contain the python code and the 
pipelines that run it respecitvely.

Pipelines
---------

* **visitQualityCore.yaml**
| The core plots for analysing the quality of the visit level data. The core pipeline is run as standard as part of the regular reprocessing. The most helpful plots go into this pipeline.

* **visitQualityExtended.yaml**
| An extended pipeline of plots and metrics to study the visit level data. The extended pipeline is run when a problem comes up or we need to look deeper into the data quality. Useful plots that address specific issues but don't need to be run on a regular basis get added to this pipeline.

* **coaddQualityCore.yaml**
| The core plots for analysing the quality of the coadd level data.

* **coaddQualityExtended.yaml**
| The extended plots for analysing the coadd data.

* **matchedVisitQualityCore.yaml**
| Plots and metrics that assess the repeatability of sources per tract by matching them between visits.

* **apCcdVisitQualityCore.yaml**
| The core plots to assess the quality of the ccd visit dataset.

python/lsst/analysis/tools
--------------------------

* **actions**
| This contains the actions that plot things and calculate things.
| Check here before adding new actions to avoid duplication.

    **scalar**

    Contains a lot of useful actions that return scalar values.
    E.g. median or sigma MAD.

    **vector**

    These actions run on vectors and return vectors.
    E.g. the S/N selector which returns an array of bools.

    **keyedData**

    These actions are base classes for other actions. You 
    shouldn't need to add stuff here. Use the scalar or 
    vector actions.

    **plots**

    The plotting code lives in here. You shouldn't need to touch 
    this unless you have to add a new plot type. Try to use one of 
    the existing ones first rather than duplicating things.

* **analysisMetrics**
| Metric classes go in here. One off metrics and very simple metrics go into analysisMetrics.py. Sets of metrics go into their own file, i.e. psfResidualMetrics.py

* **analysisParts**
| Shared code between plots and metric goes in here. Try to have as much of this as possible so that nothing changes between the plots and their associated metrics.
| I.e. shapeSizeFractionalDiff.py.

* **analysisPlots**
| Plotting classes go in here. One off plots and very simple plots go into analysisPlots.py Sets of plots go into their own file, i.e. skyObject.py.

* **contexts**
| Generic settings to be applied in a given circumstance. For example overrides that are specific to a coadd or visit level plot/metric such as the default flag selector which is different between coadd and visit analysis.

* **tasks**
| Each different dataset type requires its own task to handle the reading of the inputs.
| For example: objectTableTractAnalysis.py which handles the reading in of object tables.


A Simple Plotting Example
=========================
The first example we are going to look at is a very simple one and then we can build 
up from there. We're going to start by adapting an existing plot to our needs, we'll use a 
sky plot to show the on sky distribution of the values of a column in the table.

We use ‘actions’ to tell the code what to plot on the axis, these can be defined by anyone 
but standard ones exist already. This example will showcase some of these standard ones and 
then we’ll look more into how to define them. One of the great things about actions is that 
they allow us to only read in the columns we need from large tables.

Each plot or metric is its own class, each one has a prep, process and produce section. 
The prep section applies things like flag cuts and signal to noise cuts to the data. 
The process section builds the data required for the plot, for example if the plot 
is of a magnitude difference against a magnitude then the actions defined in the 
process section will identify which flux column needs to be read in and turned into a magnitude. 
Then another will take the fluxes needed, turn them into magnitudes and then calculate their 
difference. The produce section takes the prepared and pre calculated data and plots it on 
the graph. The plot options, such as axis labels, are set in this section.

When naming new classes it is recommended to have the word Plot in the name and that the name of the classes
matches the one that is used in the pipeline (detailed later). This name can be further expanded to include
the plot type as well.

.. code-block:: python

   class newPlot(AnalysisPlot):
       def setDefaults(self):
           super().setDefaults()
           self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
           self.prep.selectors.flagSelector.bands = ["{band}"]

           self.prep.selectors.snSelector = SnSelector()
           self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
           self.prep.selectors.snSelector.threshold = 300

           self.prep.selectors.starSelector = StarSelector()
           self.prep.selectors.starSelector.vectorKey = "{band}_extendedness"

           self.process.buildActions.xStars = LoadVector()
           self.process.buildActions.xStars.vectorKey = "coord_ra"
           self.process.buildActions.yStars = LoadVector()
           self.process.buildActions.yStars.vectorKey = "coord_dec"

           self.process.buildActions.starStatMask = SnSelector()
           self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

           self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
           self.process.buildActions.zStars.magDiff.col1 = "{band}_ap12Flux"
           self.process.buildActions.zStars.magDiff.col2 = "{band}_psfFlux"

           self.produce = SkyPlot()
           self.produce.plotTypes = ["stars"]
           self.produce.plotName = "ap12-psf_{band}"
           self.produce.xAxisLabel = "R.A. (degrees)"
           self.produce.yAxisLabel = "Dec. (degrees)"
           self.produce.zAxisLabel = "Ap 12 - PSF [mag]"
           self.produce.plotOutlines = False

Let's look at what the bits do in more detail.

.. code-block:: python

           self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
           self.prep.selectors.flagSelector.bands = ["{band}"]

The flag selector option lets us apply selectors based on flags to cut the data down. Multiple can be applied
at once and any flag that is in the input can be used. However pre built selectors already exist for the
common and recommended flag combinations.

CoaddPlotFlagSelector - this is the standard set of flags for coadd plots. The “{band}” syntax means it gets applied in the band the plot is being made in.

.. code-block:: python

           self.prep.selectors.snSelector = SnSelector()
           self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
           self.prep.selectors.snSelector.threshold = 300

SnSelector - this is the standard way of cutting the data down on S/N, you can set the flux type that is used to calculate the ratio and the threshold which the data must be above to be kept.

.. code-block:: python

           self.prep.selectors.starSelector = StarSelector()
           self.prep.selectors.starSelector.vectorKey = "{band}_extendedness"

The starSelector option is for defining a selector which picks out the specific type of object that you want
to look at. You can define this anyway you want but there are pre defined ones that can be used to choose
stars or galaxies. You can also plot both at the same time, either separately or as one dataset but the
different dynamic ranges they often cover can make the resulting plot sub optimal.

starSelector - this is the standard selector for stars. It uses the extendedness column, though any column can
be specified, the threshold in starSelector is defined for the extendedness column.

.. code-block:: python

           self.process.buildActions.xStars = LoadVector()
           self.process.buildActions.xStars.vectorKey = "coord_ra"
           self.process.buildActions.yStars = LoadVector()
           self.process.buildActions.yStars.vectorKey = "coord_dec"

This section, the xStars and yStars options, sets what is plotted on each axis. In this case it is just the
column, post selectors applied, that is directly plotted. To do this the LoadVector action is used, it just
takes a vectorKey which in this case is the column name. However this can be any action, common actions are
already defined but you can define whatever you need and use it here.

.. code-block:: python

           self.process.buildActions.starStatMask = SnSelector()
           self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

The sky plot prints some statistics on the plot, the mask that selects the points to use for these stats is
defined by the starStatMask option. In this case it uses a PSF flux based S/N selector.

.. code-block:: python

           self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
           self.process.buildActions.zStars.magDiff.col1 = "{band}_ap12Flux"
           self.process.buildActions.zStars.magDiff.col2 = "{band}_psfFlux"

The points on the sky plot are color coded by the value defined in the zStars action. Here we have gone for
the ExtinctionCorrectedMagDiff, which calculates the magnitude from each of the columns specified as col1 and
col2 and then applies extinction corrections and subtracts them. If there is no extinction corrections for the
data then it defaults to a straight difference between them.

.. code-block:: python

           self.produce = SkyPlot()
           self.produce.plotTypes = ["stars"]
           self.produce.plotName = "ap12-psf_{band}"
           self.produce.xAxisLabel = "R.A. (degrees)"
           self.produce.yAxisLabel = "Dec. (degrees)"
           self.produce.zAxisLabel = "Ap 12 - PSF [mag]"
           self.produce.plotOutlines = False

This final section declares the plot type and adds labels and things. We declare that we want to make a sky
plot, that plots only objects of type star. Next we give the plot a name that is informative for later
identification and add axis labels. The final option specifies if we want patch outlines plotted. The plot 


This new class then needs to be added to a file in analysisPlots, one off and simple plots go into the
analysisPlots file directly and the others are filed by category. For example all sky object related plots are
in the skyObjects.py file.

Once we have added the class to the relevant file we can now run it from the command line. To do this we need
to add the class to a pipeline.

.. code-block:: yaml

   description: |
     An example pipeline to run our new plot
   tasks:
     testNewPlot:
   class: lsst.analysis.tools.tasks.ObjectTableTractAnalysisTask
   config:
     connections.outputName: testNewPlot
     plots.newPlot: newPlot
   python: |
     from lsst.analysis.tools.analysisPlots import *

The class line assumes that we want to run the plot on an objectTable_tract. Each different dataset type has
its own assocaited task. Many tasks already exist for different dataset types but depending on what you want
to look at you might need to make your own.

Once we have the pipeline we can run it, the same as we would run other pipetasks.

.. code-block:: bash

   pipetask run -p pipelines/myNewPipeline.yaml
   -b /sdf/group/rubin/repo/main/butler.yaml
   -i HSC/runs/RC2/w_2022_28/DM-35609
   -o u/sr525/newPlotTest
   --register-dataset-types --prune-replaced=purge --replace-run

Let's look at each of the parts that go into the command.

.. code-block:: bash

   pipetask run -p pipelines/myNewPipeline.yaml

-p is the pipeline file, the location is relative to the directory that the command is run from.

.. code-block:: bash

   -b /sdf/group/rubin/repo/main/butler.yaml

-b is the location of the butler for the data that you want to process. This example is using the HSC data at the USDF.

.. code-block:: bash

   -i HSC/runs/RC2/w_2022_28/DM-35609

-i is the input collection to plot from, here we are using one of the weekly reprocessing runs of the RC2 data. This path is relative to the one given for the butler.yaml file in the -b option.

.. code-block:: bash

   -o u/sr525/newPlotTest

-o is the output collection that you want the plots to go into. The standard way of organising things is to put them into u/your-user-name.

.. code-block:: bash

   --register-dataset-types --prune-replaced=purge --replace-run

The other options are sometimes necessary when running the pipeline. --register-dataset-types is needed when you have a dataset type that hasn't been made before and needs to be added. --prune-replaced=purge and --replace-run are useful if you are running the same thing multiple times into the same output, for example when debugging. They replace the previous versions of the plot and just keep the most recent version.

If you don't want to include all of the data in the input collection then you need to specify a data id which
is done with the -d option.

.. code-block:: bash

   -d "instrument='HSC' AND (band='g' or band='r' or band='i' or band='z' or band='y') AND skymap='hsc_rings_v1'
   AND tract=9813 AND patch=68"

This example data id tells the processing that the instrument being used is HSC, that we want to make the plot
in the g, r, i, z and y bands, that the skymap used is the hsc_rings_v1 map, that the tract is 9813 and that
we only want to process data from patch 68 rather than all the data.


Adding an Action
================

Actions go in one of the sub folders of the actions directory depending on what type they are, this is covered in the package layout section. Before you add a new action check if it is already included before adding a duplicate. Sometimes it will probably be better to generalise an exisiting action rather than making a new one that is very similar to something that already exists. If the new action is long or specific to a given circumatance then add it to a new file, for example the ellipticity actions in ``python/lsst/analysis/tools/actions/vector/ellipticity.py``.

Let's look at some examples of actions. The first one is a scalar action.

.. code-block:: python

   class MedianAction(ScalarAction):
       vectorKey = Field[str]("Key of Vector to median.")

       def getInputSchema(self) -> KeyedDataSchema:
           return ((self.vectorKey, Vector),)

       def __call__(self, data: KeyedData, **kwargs) -> Scalar:
           mask = self.getMask(**kwargs)
           return cast(Scalar, float(np.nanmedian(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))

Let's go through what each bit of the action does.

.. code-block:: python

       vectorKey = Field[str]("Key of Vector to median.")

This is a config option, when you use the action you declare the column name using this field. This is
consistent across all actions.

.. code-block:: python

       def getInputSchema(self) -> KeyedDataSchema:
           return ((self.vectorKey, Vector),)

Every action needs a getInputSchema, this is what it uses to know which columns to read in from the table.
This means that only the needed columns can be read in allowing large tables to be accessed without memory
issues. This is one of the bonus benefits of using the ```analysis_tools``` framework.

.. code-block:: python

        def __call__(self, data: KeyedData, **kwargs) -> Scalar:
            mask = self.getMask(**kwargs)
            return cast(Scalar, float(np.nanmedian(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))

This actually does the work. It uses a mask, if it is given, and then takes the nan median of the relevant column from the data. The various calls to cast and type declarations are because it is made to work on very generic input data, any sort of keyed data type. Also we’ve got to keep typing happy otherwise we can’t merge to main.

Next we have an example of a vector action, these take vectors and return vectors.

.. code-block:: python

   class SubtractVector(VectorAction):
   """Calculate (A-B)"""

       actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
       actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

       def getInputSchema(self) -> KeyedDataSchema:
           yield from self.actionA.getInputSchema()  # type: ignore
           yield from self.actionB.getInputSchema()  # type: ignore

       def __call__(self, data: KeyedData, **kwargs) -> Vector:
           vecA = self.actionA(data, **kwargs)  # type: ignore
           vecB = self.actionB(data, **kwargs)  # type: ignore

           return vecA - vecB

Vector actions are similar to scalar actions but we will break this one down and look at the components.

.. code-block:: python

       actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
       actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

These lines are the config options, here they are the actions which give you the two values to subtract. These actions can be the loadVector action which just reads in a column without changing it in anyway.

.. code-block:: python

       def getInputSchema(self) -> KeyedDataSchema:
           yield from self.actionA.getInputSchema()  # type: ignore
           yield from self.actionB.getInputSchema()  # type: ignore

Here we get the column names from each of the actions being used, you can nest actions as deep as you want.

.. code-block:: python

       def __call__(self, data: KeyedData, **kwargs) -> Vector:
           vecA = self.actionA(data, **kwargs)  # type: ignore
           vecB = self.actionB(data, **kwargs)  # type: ignore

           return vecA - vecB

This section does the work and calculates the two actions and then subtracts them, returning the results.

These are two very simple examples of actions and how they can be used. They can be as complicated or as
simple as you want and can be composed of multiple other actions allowing common segments to be their own
actions and then reused.

Adding a Plot Type
==================

Hopefully there will be very few instances where you will need to add a new plot type and if you do please
check open ticket branches to make sure that you are not duplicating someone else's work. Try to use already
existant plot types so that we don't end up with lots of very similar plot types. Hopefully you won't really
need to touch the ploting code and can just define new classes and actions.

If you add a new plot then please make sure that you include enough providence information on the plot. There
should be enough information that anyone can recreate the plot and access the full dataset for further
investigation. See the other plots for more information on how to do this.

Current Plot Types
==================

histPlot.py -> Histograms
-------------------------
Histogram plots are constructed by the HistPlot class. These plots consist of a user-defined number of panels, with each panel consisting of a user-defined number of histograms. Panel layout will be automatically determined. A marginal table will also provide a series of summary statistics for each histogram: the total number of data points (N_data); the median of these data (Med); and the one-sigma median absolute deviation (σ_MAD). The median value will be plotted as a dashed vertical line atop the histogram, for reference.

colorColorFitPlot.py -> Stellar Locus Plots
-------------------------------------------
colorColorFitPlot is used to make the stellar locus plots. It plots a color-color diagram and then adds the
fit line from the given parameters. It also shows two panels and some statistics that detail the goodness of
the fit.

scatterPlotWithTwoHists.py -> Scatter Plots with Side Histograms
----------------------------------------------------------------
Plots a main scatter plot of the given data and then collapses the data onto a histogram for each axis. Also
calculates some statistics and adds them to the plot.

skyPlot.py -> On Sky Distribution of a Parameter
------------------------------------------------
The sky plot shows the on sky distribution of the specified value. Theoretically it can be used to plot any
x against any y with any z values to color code the points but so far the only uses have been ra against dec.

Need Help?
==========

If you get stuck with ``analysis_tools`` then feel free to reach out to the ``#rubinobs-analysis-tools``
channel on slack and hopefully someone will help you!


.. _lsst.analysis.tools-contributing:

Contributing
============


``lsst.analysis.tools`` is developed at https://github.com/lsst/analysis_tools.
You can find Jira issues for this module under the `analysis_tools <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20analysis_tools>`_ component.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 1

.. .. _lsst.analysis.tools-scripts:

.. Script reference
.. ================

.. .. TODO: Add an item to this toctree for each script reference topic in the scripts subdirectory.

.. .. toctree::
..    :maxdepth: 1

.. .. _lsst.analysis.tools-pyapi:

Python API reference
====================

.. NOTE: Skip the type definitions that cause the pipelines docs build to fail.
.. automodapi:: lsst.analysis.tools
   :no-main-docstr:
   :no-inheritance-diagram:
   :include-all-objects:
   :skip: Vector
   :skip: KeyedData
   :skip: KeyedDataSchema
   :skip: KeyedDataTypes
