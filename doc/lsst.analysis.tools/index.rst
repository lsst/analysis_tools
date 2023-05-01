.. py:currentmodule:: lsst.analysis.tools

.. _lsst.analysis.tools:

###################
lsst.analysis.tools
###################

.. Paragraph that describes what this Python module does and links to related modules and frameworks.

``analysis_tools`` is the plotting and metric framework that is used to perform QA on the pipeline products.
It is a very powerful way to explore and interact with the pipeline outputs.

.. .. _lsst.analysis.tools-using:

Using lsst.analysis.tools
=========================

For a tutorial on working with
``analysis_tools`` please see the :ref:`getting started guide <analysis-tools-getting-started>`.

.. toctree linking to topics related to using the module's APIs.

.. toctree::
   :glob:
   :maxdepth: 1

   getting-started
   action-types
   plot-types

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
.. automodapi:: lsst.analysis.tools.actions.plot
   :no-inheritance-diagram:
.. automodapi:: lsst.analysis.tools.actions.vector
   :no-inheritance-diagram:
.. automodapi:: lsst.analysis.tools.actions.scalar
   :no-inheritance-diagram:
