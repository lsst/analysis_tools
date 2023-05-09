.. py:currentmodule:: lsst.analysis.tools

.. _lsst.analysis.tools:

===================
lsst.analysis.tools
===================

The Analysis Tools package is designed to assist in the creation of quality assurance (QA) plots and metrics from the outputs of a data reduction pipeline.
The intention is that users have the flexibility to construct complex data analysis tasks from a set of simple building blocks.
In this sense, a series of consistent, repeatable and high-quality plots and metrics can be generated for any given dataset.

.. _lsst.analysis.tools-using:

Using Analysis Tools
====================

For a tutorial on working with
``analysis_tools`` please see the :ref:`getting started guide <analysis-tools-getting-started>`.

.. toctree::
   :glob:
   :maxdepth: 1

   getting-started
   action-types
   plot-types

.. _lsst.analysis.tools-help:

Need Help?
==========

If you have any questions regarding ``analysis_tools`` it is recommended that you post your question on `The Community Forum <https://community.lsst.org/>`_.

.. _lsst.analysis_tools-pyapi:

Python API Reference
====================

.. automodapi:: lsst.analysis.tools
   :include-all-objects:
.. automodapi:: lsst.analysis.tools.actions.plot
.. automodapi:: lsst.analysis.tools.actions.vector
.. automodapi:: lsst.analysis.tools.actions.scalar

.. _lsst.analysis_tools-contributing:

Contributing
============

The ``lsst.analysis.tools`` package is developed at
`github.com/lsst/analysis_tools <https://github.com/lsst/analysis_tools>`_.

Jira issues relating to this package can be found using the
`analysis_tools <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20analysis_tools>`_
component.
