#!/usr/bin/env python

# Installation script for diffpy.srfit

"""diffpy.srfit - framework for setting up complex modeling refinements.

Packages:   diffpy.srfit
Scripts:    (none yet)
"""

import os
from setuptools import setup, find_packages

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
versioncfgfile = 'diffpy/srfit/version.cfg'

def gitinfo():
    from subprocess import Popen, PIPE
    proc = Popen(['git', 'describe'], stdout=PIPE)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %ai'], stdout=PIPE)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '-'.join(desc.strip().split('-')[:2])
    rv['commit'], rv['date'] = glog.strip().split(None, 1)
    return rv


def getversioncfg():
    import os
    from ConfigParser import SafeConfigParser
    cp = SafeConfigParser()
    cp.read(versioncfgfile)
    if not os.path.isdir('.git'):  return cp
    d = cp.defaults()
    g = gitinfo()
    if g['commit'] != d.get('commit'):
        cp.set('DEFAULT', 'version', g['version'])
        cp.set('DEFAULT', 'commit', g['commit'])
        cp.set('DEFAULT', 'date', g['date'])
        cp.write(open(versioncfgfile, 'w'))
    return cp

cp = getversioncfg()

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
