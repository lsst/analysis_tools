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
If you have access to slack, project members and in kind contributers should have and get stuck with ``analysis_tools`` then feel free to reach out to the ``#rubinobs-analysis-tools``
channel on slack and hopefully someone will help you!

Another way to get help is to write a post on `The Community Forum <https://community.lsst.org/>`_.

A first place for more information is the :ref:`FAQs <analysis-tools-faqs>` page which
contains some helpful hints. More will be added to this page over time, if you find something that you think
should be added here then please do!

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

<<<<<<< HEAD
Jira issues relating to this package can be found using the
`analysis_tools <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20analysis_tools>`_
component.
=======
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
.. automodapi:: lsst.analysis.tools.actions.plot
   :no-inheritance-diagram:
.. automodapi:: lsst.analysis.tools.actions.vector
   :no-inheritance-diagram:
.. automodapi:: lsst.analysis.tools.actions.scalar
   :no-inheritance-diagram:
>>>>>>> 2de19a1 (Add more pages to the docs)
