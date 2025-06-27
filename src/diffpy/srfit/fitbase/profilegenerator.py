#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""The ProfileGenerator class for generating a profile.

ProfileGenerators encapsulate the evaluation and required Parameters and
ParameterSets of a profile calculator.  The ProfileGenerator class can
be associated with a FitContribution to help calculate a profile.

To define a ProfileGenerator, one must implement the required Parameters
and ParameterSets as well as overload the __call__ method with the
calculation. A very simple example is > class
Gaussian(ProfileGenerator): > >    def __init__(self): >        #
Initialize and give this a name >        ProfileGenerator.__init__(self,
"g") >        # Add amplitude, center and width parameters >
self.newParameter("amp", 0) >        self.newParameter("center", 0) >
self.newParameter("width", 0) > >    def __call__(self, x): >        a =
self.amp.getValue() >        x0 = self.center.getValue() >        w =
self.width.getValue() >        return a * exp(-0.5*((x-x0)/w)**2)

More examples can be found in the example directory of the
documentation.
"""

__all__ = ["ProfileGenerator"]


from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.parameterset import ParameterSet


class ProfileGenerator(Operator, ParameterSet):
    """Base class for profile generators.

    A ProfileGenerator organizes Parameters and has a __call__ method that can
    generate a profile. ProfileGenerator is also an Operator
    (diffpy.srfit.equation.literals.operators), so it can be used directly in
    an evaluation network.

    Attributes
    name            --  A name for this organizer.
    profile         --  A Profile instance that contains the calculation range
                        and will contain the generated profile.
    meta            --  A dictionary of metadata needed by the generator.
    eq              --  The Equation object used to wrap this ProfileGenerator.
                        This is set when the ProfileGenerator is added to a
                        FitContribution.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A set of constrained Parameters. Constraints can be
                        added using the 'constrain' methods.
    _parameters     --  A managed OrderedDict of contained Parameters.
    _parsets        --  A managed dictionary of ParameterSets.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    Operator Attributes
    args    --  List of Literal arguments, set with 'addLiteral'
    name    --  A name for this operator. e.g. "add" or "sin"
    nin     --  Number of inputs (<1 means this is variable)
    nout    --  Number of outputs
    operation   --  Function that performs the operation. e.g. numpy.add. In
                this case, operation is an instance method.
    symbol  --  The symbolic representation. e.g. "+" or "sin"
    _value  --  The value of the Operator.
    value   --  Property for 'getValue'.

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.
    """

    # define abstract attributes from the Operator base.
    nin = 0
    nout = 1

    def __init__(self, name):
        """Initialize the attributes."""
        Operator.__init__(self)
        ParameterSet.__init__(self, name)
        self.profile = None
        self.meta = {}
        return

    @property
    def symbol(self):
        return self.name

    # Overload me!

    def __call__(self, x):
        """Evaluate the profile.

        This method must be overloaded.

        This method only takes the independent variables to calculate
        over.
        """
        return x

    # No need to overload anything below here

    def operation(self):
        """Evaluate the profile.

        Return the result of __call__(profile.x).
        """
        y = self.__call__(self.profile.x)
        return y

    def setProfile(self, profile):
        """Assign the profile.

        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.
        """
        if self.profile is not None:
            self.profile.removeObserver(self._flush)

        self.profile = profile
        self.profile.addObserver(self._flush)
        self._flush(other=(self,))

        # Merge the profiles metadata with our own
        self.meta.update(self.profile.meta)
        self.processMetaData()
        return

    def processMetaData(self):
        """Process the metadata.

        This can be used to configure a ProfileGenerator upon a change
        in the metadata. This method gets called whenever the Profile is
        set.
        """
        return

    def _validate(self):
        """Validate my state.

        This performs profile validations. This performs ParameterSet
        validations. This does not validate the operation, since this
        could be costly. The operation should be validated with a
        containing equation.

        Raises SrFitError if validation fails.
        """
        if self.profile is None:
            raise SrFitError("profile is None")
        self.profile._validate()
        ParameterSet._validate(self)
        return


# End class ProfileGenerator
