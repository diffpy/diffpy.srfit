.. _developers-guide-index:

####################################################
diffpy.srfit documentation
####################################################

diffpy.srfit - configurable code for solving atomic structures.

| Software version |release|.
| Last updated |today|.

The diffpy.srfit package provides the framework for building a global optimizer
on the fly from components such as function calculators (that calculate
different data spectra), regression algorithms and structure models.  The
software is capable of co-refinement using multiple information sources or
models. It provides a uniform interface for various regression algorithms. The
target function being optimized can be specified by the user according to the
data available.

Within the diffpy.srfit framework, any parameter used in describing the
structure of a material can be passed as a refinable variable to the global
optimizer.  Once parameters are declared as variables they can easily be turned
"on" or "off", i.e. fixed or allowed to vary. Additionally, variables may be
constrained to obey mathematical relationships with other parameters or
variables used in the structural model. Restraints can be applied to
variables, which adds a penalty to the refinement process commensurate with the
deviation from the known value or range. The cost function can also be
customized by the user. If the refinement contains multiple models, each model
can have its own cost function which will be properly weighted and combined to
obtain the total cost function. Additionally, diffpy.srfit is designed to be
extensible, allowing the user to integrate external calculators to perform
co-refinements with other techniques.

========================================
Authors
========================================

diffpy.srfit is developed by members of the Billinge Group at
Columbia University and at Brookhaven National Laboratory including
Christopher L. Farrow, Pavol Juhás, Simon J.L. Billinge.

The source code in *observable.py* was derived from the 1.0 version
of the Caltech "Pyre" project.

For a detailed list of contributors see
https://github.com/diffpy/diffpy.srfit/graphs/contributors.

======================================
Installation
======================================

See the `README <https://github.com/diffpy/diffpy.srfit#requirements>`_
file included with the distribution.

======================================
Where next?
======================================

.. toctree::
   :maxdepth: 2

   examples.rst
   extending.rst

======================================
Table of contents
======================================

.. toctree::
   :titlesonly:

   license
   api/diffpy.srfit.rst

.. faq.rst

* :ref:`genindex`
* :ref:`search`
