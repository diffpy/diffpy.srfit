#!/usr/bin/env python

# Installation script for diffpy.srfit

"""diffpy.srfit - framework for setting up complex modeling refinements.

Packages:   diffpy.srfit
Scripts:    (none yet)
"""

import os
from setuptools import setup, find_packages

def gitversion():
    from subprocess import Popen, PIPE
    proc = Popen(['git', 'describe'], stdout=PIPE)
    desc = proc.stdout.read().strip()
    proc = Popen(['git', 'log', '-1', '--format=%ai'], stdout=PIPE)
    isodate = proc.stdout.read()
    date = isodate.split()[0].replace('-', '')
    rv = desc + '-' + date
    return rv


def getsetupcfg():
    cfgfile = 'setup.cfg'
    from ConfigParser import SafeConfigParser
    cp = SafeConfigParser()
    cp.read(cfgfile)
    if not os.path.isdir('.git'):  return cp
    d = cp.defaults()
    vcfg = d.get('version', '')
    vgit = gitversion()
    if vgit != vcfg:
        cp.set('DEFAULT', 'version', vgit)
        cp.write(open(cfgfile, 'w'))
    return cp

cp = getsetupcfg()

# define distribution
dist = setup(
        name = "diffpy.srfit",
        version = cp.get('DEFAULT', 'version'),
        namespace_packages = ['diffpy'],
        packages = find_packages(exclude=['tests']),
        test_suite = 'diffpy.srfit.tests',
        include_package_data = True,
        entry_points = {},
        install_requires = [
            'diffpy.Structure>=1.0-r5333-20100518',
            'pyobjcryst>=1.0b1.dev-r5681-20100816',
            'diffpy.srreal>=0.2a1.dev-r6037-20101130',
            'periodictable>=1.0b1.dev-r5681-20100816',
            'numpy>=1.0',
            'scipy>=0.7.0',
            ],
        dependency_links = [
            # REMOVE dev.danse.us for a public release.
            'http://dev.danse.us/packages/',
            "http://www.diffpy.org/packages/",
        ],

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
