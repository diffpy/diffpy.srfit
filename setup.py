#!/usr/bin/env python

# Installation script for diffpy.Structure

"""diffpy.srfit - prototype for new PDF calculator and assortment
of real space utilities.

Packages:   diffpy.srfit
Scripts:    (none yet)
"""

from setuptools import setup, find_packages
from setuptools import Extension
import fix_setuptools_chmod
import glob

clickermodule = Extension('diffpy.srfit.util._clicker', 
        glob.glob("extensions/util/*.cpp"),
        libraries = ["boost_python"],
        )

# define distribution
dist = setup(
        name = "diffpy.srfit",
        version = "1.0a4",
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        entry_points = {},
        install_requires = [
            'diffpy.Structure',
            'pyobjcryst>=0.1a1.dev-r3762',
            'periodictable>=0.9.dev-r197-20090419',
            ],
        dependency_links = [
            # REMOVE dev.danse.us for a public release.
            'http://dev.danse.us/packages/',
            "http://www.diffpy.org/packages/",
        ],

        # The extensions
        ext_modules = [clickermodule],

        # Stuff for building extensions
        zip_safe = False,

        author = "Simon J.L. Billinge",
        author_email = "sb2896@columbia.edu",
        maintainer = 'Christopher L. Farrow',
        maintainer_email = 'clf2121@columbia.edu',
        description = "SrFit - Structure refinement from diffraction data",
        license = "BSD",
        url = "http://www.diffpy.org/",
        keywords = "complex modeling calculator utilities",
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 2.6',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
