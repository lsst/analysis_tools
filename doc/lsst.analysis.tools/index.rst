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
   faqs
   detailed-package-layout

.. _lsst.analysis.tools-help:

Need Help?
==========

The first place to check for more information is the :ref:`FAQs <analysis-tools-faqs>` page which contains some helpful hints.
More will be added to this page over time; if you find something that you think should be added here then please do!

Project members and in-kind contributors may ask for further help via the ``#rubinobs-analysis-tools`` channel on the LSSTC Slack workspace.
All users of the analysis tools package, including community members, are also welcome to post any questions you might have on the `Community Forum <https://community.lsst.org/>`_.

.. _lsst.analysis_tools-contributing:

Contributing
============

The ``lsst.analysis.tools`` package is developed at `github.com/lsst/analysis_tools <https://github.com/lsst/analysis_tools>`_.

Jira issues relating to this package can be found using the `analysis_tools <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20analysis_tools>`_ component.

.. _lsst.analysis_tools-scripts:

Script reference
================

.. toctree::
   :maxdepth: 1

   scripts/verify_to_sasquatch.py
   scripts/build-gather-resource-usage-qg

.. _lsst.analysis_tools-pyapi:

Python API Reference
====================

.. automodapi:: lsst.analysis.tools
   :include-all-objects:
.. automodapi:: lsst.analysis.tools.actions.plot
.. automodapi:: lsst.analysis.tools.actions.scalar
.. automodapi:: lsst.analysis.tools.actions.vector
.. automodapi:: lsst.analysis.tools.actions.keyedData
.. automodapi:: lsst.analysis.tools.tasks
