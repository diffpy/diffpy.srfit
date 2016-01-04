#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2015 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################


"""\
Universal import functions for volatile SasView/SansViews API-s.
"""

def sasimport(modname):
    """Import specified module from the SasView sas package.

    modname  -- absolute module name contained in the sas package.

    When specified import does not work directly, try older API-s and raise
    DeprecationWarning.  Raise ImportError if nothing works.
    """
    if not modname.startswith('sas.'):
        emsg = 'Module name must start with "sas."'
        raise ValueError(emsg)
    mobj = None
    # standard import
    try:
        exec('import %s as mobj' % modname)
    except ImportError:
        pass
    else:
        return mobj
    # revert to the old sans namespace, sas --> sans
    modsans = 'sans' + modname[3:]
    import warnings
    wfmt = ("Using obsolete package %r instead of %r.  Please install "
            "SasView 3.1 or the srfit-sasview package from Anaconda.")
    wmsg = wfmt % (modsans, modname)
    try:
        exec('import %s as mobj' % modsans)
        warnings.warn(wmsg, DeprecationWarning)
    except ImportError:
        pass
    else:
        return mobj
    # finally check the oldest DataLoader API for sas.dataloader
    if modname.startswith('sas.dataloader'):
        modloader = 'DataLoader' + modname[14:]
        wmsg = wfmt % (modloader, modname)
        try:
            exec('import %s as mobj' % modloader)
            warnings.warn(wmsg, DeprecationWarning)
        except ImportError:
            pass
        else:
            return mobj
    # Obsolete API-s failed here.  Import again and let it raise ImportError.
    exec('import %s as mobj' % modname)
    raise AssertionError("The import above was supposed to fail.")
