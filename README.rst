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
constrained to obey mathematical relationsips with other parameters or
variables used in the structural model. Restraints can alo be applied to
variables, which adds a penalty to the refinment process commensurate with the
deviation from the known value or range. The cost function can also be
customized by the user. If the refinement contains multiple models, each model
can have its own cost function which will be properly weighted and combined to
obtain the total cost function. Additionaly, diffpy.srfit is designed to be
extensible, allowing the user to integrate external calculators to perform
corefinements with other techniques.

For more information about the diffpy.srfit library, see the users manual at
http://diffpy.github.io/diffpy.srfit.

REQUIREMENTS
------------------------------------------------------------------------

The diffpy.srfit package requires Python 2.6 or 2.7 and the following software:

* ``setuptools`` - software distribution tools for Python
* ``NumPy`` - numerical mathematics and fast array operations for Python
* ``SciPy`` - scientific libraries for Python
* ``matplotlib`` - python plotting library

The functions in the diffpy.srfit.structure module require:

* ``diffpy.Structure`` - crystal structure container and parsers
* ``pyobjcryst`` - python bindings to ObjCryst++
* ``diffpy.srreal`` - python library for PDF calculation

The functions in the diffpy.srfit.pdf module require:

* ``diffpy.srreal`` - python library for PDF calculation

The functions in the diffpy.srfit.sas module require:

* ``sans.pr`` - calculation of P(R) in small-angle scattering

On Ubuntu Linux, the required software can easily be installed using
the system package manager::

    sudo aptitude install python-setuptools python-numpy python-scipy python-matplotlib

For Mac OS X systems with the MacPorts package manager, the required
software can be installed with::

    sudo port install python27 py27-setuptools py27-numpy py27-scipy py27-matplotlib

When installing for MacPorts, make sure the MacPorts bin directory is the first
in the system PATH and that python27 is selected as the default Python version
in MacPorts::

    sudo port select --set python python27

For other required packages see their respective web pages for installation
instructions.


INSTALLATION
------------------------------------------------------------------------

The easiest option is to use the latest DiffPy-CMI release bundle from
http://www.diffpy.org/, which comes with diffpy.srfit and all other
dependencies included.

Or, use ``easy_install`` to download and install the latest release from
`Python Package Index <https://pypi.python.org>`_::

    sudo easy_install diffpy.srfit

If you prefer to install from sources, make sure all required software packages
are in place and then run::

    sudo python setup.py install

This installs diffpy.srfit for all users in the default system location.  If
administrator (root) access is not available, see the usage info from ``python
setup.py install --help`` for options to install to user-writable directories.
The installation integrity can be verified by changing to the HOME directory
and running::

    python -m diffpy.srfit.tests.run


DEVELOPMENT
------------------------------------------------------------------------

diffpy.srfit is an open-source software developed as a part of the DiffPy-CMI
complex modeling initiative at the Brookhaven National Laboratory.  The
diffpy.srfit sources are hosted at
https://github.com/diffpy/diffpy.srfit/,

Feel free to fork the project and contribute.  To install diffpy.srfit in a
development mode, with its sources being directly used by Python rather than
copied to a package directory, use::

    python setup.py develop --user


ACKNOWLEDGEMENT
------------------------------------------------------------------------

Part of the source code in *_abc.py* and *_ordereddict.py* was derived from
python2.7 at http://svn.python.org/; while other code *observable.py*
was derived the 1.0 version of the Caltech 'pyre' project.


CONTACTS
------------------------------------------------------------------------

For more information on diffpy.srfit please visit the project web-page

http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu.
