
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
"""Building blocks for defining a state-aware equation

Literals are the building blocks for equations. Literals are composed into
equation trees (or Literal trees) that can be evaluated or have other actions
performed on them by Visitors (in diffpy.srfit.equation.visitors). The
Literal-Visitor relationship is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern), execept that Literals contain
attributes that are used specifically by the Evaluator visitor.

The simplest Literal is the Argument, which holds the name and value of a
equation variable.  Operators are Literals that are used to compose other
Literals. Once an Operator is evaluated, it has a value that it can pass to
other operators. Hence, a network or tree of Literals can be composed that
evaluate as an equation. 

Literals have a Clicker (diffpy.srfit.equation.Clicker) member that can inform
its direct ancestor Literals about a change in its value. These ancestors
Literals inform their ancestors of the change as well, and in this way all
ancestors of an Argument know when it changes. Visitors can contain their own
clicker that can be compared to that of a Literal, and make a decision on
whether to reprocess the information on a given node of the Literal tree. 

Modules:
Argument    --  Contains the Argument class (see below). The module
                documentation contains more information about Arguments.
Generator   --  Contains the Generator class (see below). The module
                documentation contains more information about Generators.
Literal     --  Contains the Literal base class. The module documentation
                contains more information about Literals.
operators   --  Contains various Operator classes (see below). The module
                documentation contains more information about Operators.
Partition   --  Contains the Partition class (see below). The module
                documentation contains more information about Partitions.

Classes:
Argument    --  Contains the Argument class. Arguments have a value and name,
                and can be made to appear as constant. Arguments are leaf nodes
                of an equation tree.
Generator   --  A Generator is a Literal that can generate other Literals. Use
                of generators help extend the architecture to include equations
                that are not easily expressed as with other Literals.
                Generators are evaluated like leaf nodes of an equation tree,
                but can be dependent on other Literals.
Operator    --  An Operator can produce a new value based on the values of the
                Literals it contains. Operators can contain tags that tell it
                how to operate on tagged Partitions. Operators are internal
                nodes of a Literal tree.
                Operators:
                AdditionOperator        --  Encapsulates Addition
                SubtractionOperat       --  Encapsulates Subtraction
                MultiplicationOperato   --  Encapsulates Multiplication
                DivisionOperator        --  Encapsulates Division
                ExponentiationOperato   --  Encapsulates Exponentiation
                RemainderOperator       --  Encapsulates Remainder (Modulo)
                NegationOperator        --  Encapsulates Negation
                UFuncOperator           -- Wraps a numpy ufunc
                CombineOperator         -- Explicitly combines a partition
Partition   --  Partitions are Literals that explicitly decompose a value into
                logical pieces, and can recompose these into a single value.
                Partitions are leaf nodes in the Literal tree, but are
                fundamentally a collection of tagged Arguments. The tags are
                are used by Operators for conditional operation on different
                pieces of the Partition.

"""

# package version
from diffpy.srfit.version import __version__

# Import the operators

from .Argument import Argument
from .Generator import Generator
from .Partition import Partition
from .operators import Operator
from .operators import AdditionOperator
from .operators import SubtractionOperator
from .operators import MultiplicationOperator
from .operators import DivisionOperator
from .operators import ExponentiationOperator
from .operators import RemainderOperator
from .operators import NegationOperator
from .operators import UFuncOperator
from .operators import CombineOperator
