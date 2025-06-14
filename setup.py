#!/usr/bin/env python

# Installation script for diffpy.srfit
"""diffpy.srfit - framework for setting up complex modeling refinements.

Packages:   diffpy.srfit
"""

import os
import re
import sys

from setuptools import find_packages, setup

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = "3.0.0.post0"

# determine if we run with Python 3.
PY3 = sys.version_info[0] == 3

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
MYDIR = os.path.dirname(os.path.abspath(__file__))
versioncfgfile = os.path.join(MYDIR, "src/diffpy/srfit/version.cfg")
gitarchivecfgfile = os.path.join(MYDIR, ".gitarchive.cfg")


def gitinfo():
    from subprocess import PIPE, Popen

    kw = dict(stdout=PIPE, cwd=MYDIR, universal_newlines=True)
    proc = Popen(["git", "describe", "--match=v[[:digit:]]*"], **kw)
    desc = proc.stdout.read()
    proc = Popen(["git", "log", "-1", "--format=%H %ct %ci"], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv["version"] = ".post".join(desc.strip().split("-")[:2]).lstrip("v")
    rv["commit"], rv["timestamp"], rv["date"] = glog.strip().split(None, 2)
    return rv


def getversioncfg():
    if PY3:
        from configparser import RawConfigParser
    else:
        from ConfigParser import RawConfigParser
    vd0 = dict(version=FALLBACK_VERSION, commit="", date="", timestamp=0)
    # first fetch data from gitarchivecfgfile, ignore if it is unexpanded
    g = vd0.copy()
    cp0 = RawConfigParser(vd0)
    cp0.read(gitarchivecfgfile)
    if len(cp0.get("DEFAULT", "commit")) > 20:
        g = cp0.defaults()
        mx = re.search(r"\btag: v(\d[^,]*)", g.pop("refnames"))
        if mx:
            g["version"] = mx.group(1)
    # then try to obtain version data from git.
    gitdir = os.path.join(MYDIR, ".git")
    if os.path.exists(gitdir) or "GIT_DIR" in os.environ:
        try:
            g = gitinfo()
        except OSError:
            pass
    # finally, check and update the active version file
    cp = RawConfigParser()
    cp.read(versioncfgfile)
    d = cp.defaults()
    rewrite = not d or (g["commit"] and (g["version"] != d.get("version") or g["commit"] != d.get("commit")))
    if rewrite:
        cp.set("DEFAULT", "version", g["version"])
        cp.set("DEFAULT", "commit", g["commit"])
        cp.set("DEFAULT", "date", g["date"])
        cp.set("DEFAULT", "timestamp", g["timestamp"])
        with open(versioncfgfile, "w") as fp:
            cp.write(fp)
    return cp


versiondata = getversioncfg()

with open(os.path.join(MYDIR, "README.rst")) as fp:
    long_description = fp.read()

# define distribution
setup_args = dict(
    name="diffpy.srfit",
    version=versiondata.get("DEFAULT", "version"),
    packages=find_packages(os.path.join(MYDIR, "src")),
    package_dir={"": "src"},
    test_suite="diffpy.srfit.tests",
    include_package_data=True,
    install_requires=["six"],
    zip_safe=False,
    author="Simon J.L. Billinge",
    author_email="sb2896@columbia.edu",
    maintainer="Pavol Juhas",
    maintainer_email="pavol.juhas@gmail.com",
    description="SrFit - Structure refinement from diffraction data",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="BSD-style license",
    url="https://github.com/diffpy/diffpy.srfit",
    keywords="optimization constraints restraints structure refinement complex modeling",
    classifiers=[
        # List of possible values at
        # http://pypi.python.org/pypi?:action=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
    ],
)

if __name__ == "__main__":
    setup(**setup_args)

# End of file
