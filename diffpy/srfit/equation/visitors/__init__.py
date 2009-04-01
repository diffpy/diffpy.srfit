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
"""Visitors that perform on Literal trees.

Visitors are designed to traverse and extract infromation from Literal trees
(diffpy.srfit.equation.literals). The primary Visitor is the Evaluator, which
can evaluate the value of a Literal tree. The Evaluator can inspect the state
of a Literal tree and decide if it must re-evaluate values of given nodes on
subsequent visits. Other Visitors are useful for validating, printing and
extracting Argument literals from the Literal tree.

The Literal-Visitor relationship is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern), execept that Literals contain
attributes that are used specifically by the Evaluator visitor.

Modules:
ArgFinder   --  Contains the ArgFinder class (see below). The module
                documentation contains more information about ArgFinders.
Evaluator   --  Contains the Evaluator class (see below). The module
                documentation contains more information about Evaluators.
Printer     --  Contains the Printer class (see below). The module
                documentation contains more information about Printers.
Validiator  --  Contains the Validator class (see below). The module
                documentation contains more information about Validators.
Visitor     --  Contains the Visitor base class. The module documentation
                contains more information about Visitors.

Classes:
ArgFinder   --  ArgFinder extracts Arguments from a Literal tree.
Evaluator   --  Evaluates a Literal tree.
Printer     --  A not-so-pretty-printer for an Argument tree. Used primarily
                for debugging.
Validator   --  Used to validate a Literal tree.

"""

# package version
from diffpy.srfit.version import __version__

from .ArgFinder import ArgFinder
from .Printer import Printer
from .Evaluator import Evaluator
from .Validator import Validator
