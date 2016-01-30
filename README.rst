.. image:: https://travis-ci.org/diffpy/diffpy.srfit.svg?branch=master
   :target: https://travis-ci.org/diffpy/diffpy.srfit

.. image:: http://codecov.io/github/diffpy/diffpy.srfit/coverage.svg?branch=master
   :target: http://codecov.io/github/diffpy/diffpy.srfit?branch=master


diffpy.srfit
========================================================================

Configurable code for solving atomic structures.

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

For more information about the diffpy.srfit library, see the users manual at
http://diffpy.github.io/diffpy.srfit.

REQUIREMENTS
------------------------------------------------------------------------

The diffpy.srfit package requires Python 2.7 and the following software:

* ``setuptools`` - software distribution tools for Python
* ``NumPy`` - numerical mathematics and fast array operations for Python
* ``SciPy`` - scientific libraries for Python
* ``matplotlib`` - python plotting library

Recommended software:

Optimizations involving crystal structures or molecules require

* ``diffpy.Structure`` - crystal structure container and parsers,
  https://github.com/diffpy/diffpy.Structure
* ``pyobjcryst`` - Crystal and Molecule storage, rigid units, bond
  length and bond angle restraints, https://github.com/diffpy/pyobjcryst

Optimizations involving pair distribution functions PDF or bond valence
sums require

* ``diffpy.srreal`` - python library for PDF calculation,
  https://github.com/diffpy/diffpy.srreal

Optimizations involving small angle scattering or shape characteristic
functions from the diffpy.srfit.sas module require

* ``sas`` - module for calculation of P(R) in small-angle scattering
  from the SasView project, http://www.sasview.org

We recommend to use `Anaconda Python <https://www.continuum.io/downloads>`_
as it allows to install all software dependencies together with
diffpy.srfit.  For other Python distributions it is necessary to
install the required software separately.  As an example, on Ubuntu
Linux some of the required software can be installed using ::

   sudo apt-get install \
      python-setuptools python-numpy python-scipy python-matplotlib

For other required packages see their respective web pages for installation
instructions.


INSTALLATION
------------------------------------------------------------------------

The preferred method is to use Anaconda Python and install from the
"diffpy" channel of Anaconda packages ::

   conda config --add channels diffpy
   conda install diffpy.srfit

diffpy.srfit is also included in the "diffpy-cmi" collection
of packages for structure analysis ::

   conda install diffpy-cmi

Another option is to use ``easy_install`` to download and install the
latest release from `Python Package Index <https://pypi.python.org>`_ ::

   easy_install diffpy.srfit

If you prefer to install from sources, make sure all required software
packages are in place and then run ::

   python setup.py install

You may need to use ``sudo`` with system Python so the process is
allowed to put files to the system directories.  If administrator (root)
access is not available, consult the output from
``python setup.py install --help`` for options to install to a
user-writable locations.  The installation integrity can be verified by
changing to the HOME directory and running ::

   python -m diffpy.srfit.tests.run


DEVELOPMENT
------------------------------------------------------------------------

diffpy.srfit is an open-source software developed as a part of the DiffPy-CMI
complex modeling initiative at the Brookhaven National Laboratory.  The
diffpy.srfit sources are hosted at
https://github.com/diffpy/diffpy.srfit.

Feel free to fork the project and contribute.  To install diffpy.srfit
in a development mode, with its sources being directly used by Python
rather than copied to a package directory, use ::

   python setup.py develop --user


ACKNOWLEDGEMENT
------------------------------------------------------------------------

Part of the source code in *_abc.py* and *_ordereddict.py* was derived from
Python 2.7 at http://www.python.org/download/source; while other code
*observable.py* was derived from the 1.0 version of the Caltech "Pyre"
project.


CONTACTS
------------------------------------------------------------------------

For more information on diffpy.srfit please visit the project web-page

http://www.diffpy.org

or email Prof. Simon Billinge at sb2896@columbia.edu.
