#!/usr/bin/env python
##############################################################################
#
# (c) 2008-2025 The Trustees of Columbia University in the City of New York.
# All rights reserved.
# (c) 2026-present The DiffPy Team. All rights reserved.
#
# File coded by: Christopher Farrow, Pavol Juhas, Caden Myers,
# Simon J. L. Billinge, and members of the DiffPy community.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.srfit/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Definition of __version__."""

#  We do not use the other three variables, but can be added back if needed.
#  __all__ = ["__date__", "__git_commit__", "__timestamp__", "__version__"]

# obtain version information
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("diffpy.srfit")
except PackageNotFoundError:
    __version__ = "unknown"
