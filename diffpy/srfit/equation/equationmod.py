#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""The Equation class for holding and evaluating an equation.

Equation is a functor that holds a Literal tree that defines an equation. It's
__call__ method evaluates the equation at the most recent value of its
Arguments. The non-constant arguments are accessible as attributes of the
Equation instance.

"""

__all__ = ["Equation"]


from itertools import ifilter
from types import ModuleType

import numpy

from argument import Argument
from diffpy.srfit.util.ordereddict import OrderedDict

__f = lambda x : not (isinstance(x[1], ModuleType) or x[0].startswith("_"))
ns = dict( ifilter(__f, numpy.__dict__.iteritems() ) )

class Equation(object):
    """Class for holding and evaluating an Equation."""

    def __init__(self, eqstr=""):
        """Initialize.
        """
        self.eqstr = eqstr
        self.args = OrderedDict()
        self.ops = {}

        return

    def __getattr__(self, name):
        """Gives access to the Arguments as attributes."""
        arg = self.args.get(name)
        if arg is None:
            raise AttributeError("No argument named '%s' here"%name)
        return arg

    def addArg(self, arg):
        self.args[arg.name] = arg
        return

    def addOp(self, op, name):
        self.ops[name] = op
        return

    def setEqString(self, eqstr):
        """Set the root of the Literal tree.

        Raises:
        ValueError if errors are found in the Literal tree.

        """
        self.eqstr = eqstr
        return

    def __call__(self, *args, **kw):
        """Call the equation.
        
        New Argument values are acceped as arguments or keyword assignments (or
        both).  The order of accepted arguments is given by the args
        attribute.  The Equation will remember values set in this way.

        Raises
        ValueError when a passed argument cannot be found

        """

        # Process args
        if len(args) != len(self.args):
            raise ValueError("Too many arguments")
        if args:
            locals = dict(zip(self.args.keys(), args))
        else:
            locals = dict([(k, v.getValue()) for k, v in self.args.items()])

        # Process kw
        locals.update(kw)

        locals.update(self.ops)

        # Evaluate the function
        return eval(self.eqstr, ns, locals)


# version
__id__ = "$Id$"

#
# End of file
