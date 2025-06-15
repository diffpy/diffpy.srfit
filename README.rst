|Icon| |title|_
===============

.. |title| replace:: diffpy.srfit
.. _title: https://diffpy.github.io/diffpy.srfit

.. |Icon| image:: https://avatars.githubusercontent.com/diffpy
        :target: https://diffpy.github.io/diffpy.srfit
        :height: 100px

|PyPI| |Forge| |PythonVersion| |PR|

|CI| |Codecov| |Black| |Tracking|

.. |Black| image:: https://img.shields.io/badge/code_style-black-black
        :target: https://github.com/psf/black

.. |CI| image:: https://github.com/diffpy/diffpy.srfit/actions/workflows/matrix-and-codecov-on-merge-to-main.yml/badge.svg
        :target: https://github.com/diffpy/diffpy.srfit/actions/workflows/matrix-and-codecov-on-merge-to-main.yml

.. |Codecov| image:: https://codecov.io/gh/diffpy/diffpy.srfit/branch/main/graph/badge.svg
        :target: https://codecov.io/gh/diffpy/diffpy.srfit

.. |Forge| image:: https://img.shields.io/conda/vn/conda-forge/diffpy.srfit
        :target: https://anaconda.org/conda-forge/diffpy.srfit

.. |PR| image:: https://img.shields.io/badge/PR-Welcome-29ab47ff

.. |PyPI| image:: https://img.shields.io/pypi/v/diffpy.srfit
        :target: https://pypi.org/project/diffpy.srfit/

.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/diffpy.srfit
        :target: https://pypi.org/project/diffpy.srfit/

.. |Tracking| image:: https://img.shields.io/badge/issue_tracking-github-blue
        :target: https://github.com/diffpy/diffpy.srfit/issues

diffpy.srfit
============

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

For more information about the diffpy.srfit library, please consult our `online documentation <https://diffpy.github.io/diffpy.srfit>`_.

Citation
--------

If you use diffpy.srfit in a scientific publication, we would like you to cite this package as


   P. Juhás, C. L. Farrow, X. Yang, K. R. Knox and S. J. L. Billinge,
   `Complex modeling: a strategy and software program for combining
   multiple information sources to solve ill posed structure and
   nanostructure inverse problems
   <http://dx.doi.org/10.1107/S2053273315014473>`__,
   *Acta Crystallogr. A* **71**, 562-568 (2015).

Installation
------------

The preferred method is to use `Miniconda Python
<https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>`_
and install from the "conda-forge" channel of Conda packages.

To add "conda-forge" to the conda channels, run the following in a terminal. ::

        conda config --add channels conda-forge

We want to install our packages in a suitable conda environment.
The following creates and activates a new environment named ``diffpy.srfit_env`` ::

        conda create -n diffpy.srfit_env diffpy.srfit
        conda activate diffpy.srfit_env

To confirm that the installation was successful, type ::

        python -c "import diffpy.srfit; print(diffpy.srfit.__version__)"

The output should print the latest version displayed on the badges above.

This will install the minimal `diffpy.srfit` installation.  It will often be used
as along with other packages for manipulating and computing crystal structures
and so on.  We also therefore recommend installing the following:

* ``diffpy.structure`` - crystal structure container and parsers,
  https://github.com/diffpy/diffpy.structure
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

Fallback Installation
---------------------

If the above does not work, you can use ``pip`` to download and install the latest release from
`Python Package Index <https://pypi.python.org>`_.
To install using ``pip`` into your ``diffpy.srfit_env`` environment, type ::

        pip install diffpy.srfit

If you prefer to install from sources, after installing the dependencies, obtain the source archive from
`GitHub <https://github.com/diffpy/diffpy.srfit/>`_. Once installed, ``cd`` into your ``diffpy.srfit`` directory
and run the following ::

        pip install .

Getting Started
---------------

You may consult our `online documentation <https://diffpy.github.io/diffpy.srfit>`_ for tutorials and API references.

Support and Contribute
----------------------

If you see a bug or want to request a feature, please `report it as an issue <https://github.com/diffpy/diffpy.srfit/issues>`_ and/or `submit a fix as a PR <https://github.com/diffpy/diffpy.srfit/pulls>`_.

Feel free to fork the project and contribute. To install diffpy.srfit
in a development mode, with its sources being directly used by Python
rather than copied to a package directory, use the following in the root
directory ::

        pip install -e .

To ensure code quality and to prevent accidental commits into the default branch, please set up the use of our pre-commit
hooks.

1. Install pre-commit in your working environment by running ``conda install pre-commit``.

2. Initialize pre-commit (one time only) ``pre-commit install``.

Thereafter your code will be linted by black and isort and checked against flake8 before you can commit.
If it fails by black or isort, just rerun and it should pass (black and isort will modify the files so should
pass after they are modified). If the flake8 test fails please see the error messages and fix them manually before
trying to commit again.

Improvements and fixes are always appreciated.

Before contributing, please read our `Code of Conduct <https://github.com/diffpy/diffpy.srfit/blob/main/CODE_OF_CONDUCT.rst>`_.

Contact
-------

For more information on diffpy.srfit please visit the project `web-page <https://diffpy.github.io/>`_ or email Simon Billinge at sb2896@columbia.edu.

Acknowledgements
----------------

``diffpy.srfit`` is built and maintained with `scikit-package <https://scikit-package.github.io/scikit-package/>`_.

The source code in *observable.py* was derived from the 1.0 version
of the Caltech "Pyre" project.
