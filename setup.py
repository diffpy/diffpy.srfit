#!/usr/bin/env python

# Installation script for diffpy.Structure

"""diffpy.srfit - prototype for new PDF calculator and assortment
of real space utilities.

Packages:   diffpy.srfit
Scripts:    (none yet)
"""

from setuptools import setup, find_packages
import fix_setuptools_chmod

# define distribution
dist = setup(
        name = "diffpy.srfit",
        version = "1.0a2",
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        entry_points = {},
        install_requires = [],
        dependency_links = [
            # REMOVE dev.danse.us for a public release.
            'http://dev.danse.us/packages/',
            "http://www.diffpy.org/packages/",
        ],

        author = "Simon J.L. Billinge",
        author_email = "sb2896@columbia.edu",
        description = "Package for structure refinement from diffraction data",
        license = "BSD",
        url = "http://www.diffpy.org/",
        keywords = "complex modeling calculator utilities",
)

# End of file
