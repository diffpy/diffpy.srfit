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
calculation. A very simple example is

.. code:: python

    class Gaussian(ProfileGenerator):
        def __init__(self):
            # Initialize and give this a name
            ProfileGenerator.__init__(self, "g")
            # Add amplitude, center and width parameters
            self.new_parameter("amp", 0)
            self.new_parameter("center", 0)
            self.new_parameter("width", 0)
        def __call__(self, x):
            a = self.amp.get_value()
            x0 = self.center.get_value()
            w = self.width.get_value()
            return a * exp(-0.5*((x-x0)/w)**2)

More examples can be found in the example directory of the
documentation.
"""

__all__ = ["ProfileGenerator"]


from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.parameterset import ParameterSet


class ProfileGenerator(Operator, ParameterSet):
    """Base class for profile generators.

    A ProfileGenerator organizes Parameters and has a ``__call__`` method
    that can generate a profile. ProfileGenerator is also an Operator
    (diffpy.srfit.equation.literals.operators), so it can be used directly
    in an evaluation network.

    Parameters
    ----------
    name : str
        The name of the ProfileGenerator. It must be a valid attribute
        name.

    Attributes
    ----------
    name : str
        The name of the ProfileGenerator.
    profile : Profile or None
        The Profile that contains the calculation range and will contain
        the generated profile. It is None until ``set_profile`` is called.
    meta : dict
        The metadata needed by the generator. It is updated from the
        Profile metadata whenever the Profile is set or changes.
    names : list of str
        The names of the top-level Parameters (read only). See
        ``get_names``.
    values : list
        The values of the top-level Parameters (read only). See
        ``get_values``.
    args : list
        The Literal arguments of the operator, added with
        ``add_literal``.
    nin : int
        The number of inputs, which is 0. The independent variable is
        taken from ``profile``.
    nout : int
        The number of outputs, which is 1.
    symbol : str
        The symbolic representation of the operator, which is the same as
        ``name`` (read only).
    value : object
        The value of the operator (read only). See ``get_value``.
    _value : object
        The cached value of the operator, or None if it must be
        recomputed.
    _calculators : dict
        The managed dictionary of Calculators, indexed by name.
    _constraints : dict
        The dictionary of Constraints, indexed by the constrained
        Parameter. Constraints can be added using the ``add_constraint``
        method.
    _parameters : OrderedDict
        The managed OrderedDict of contained Parameters.
    _parsets : OrderedDict
        The managed OrderedDict of contained ParameterSets.
    _restraints : set
        The set of Restraints. Restraints can be added using the
        ``add_soft_bounds`` method.
    _eqfactory : diffpy.srfit.equation.builder.EquationFactory
        The factory that is used to create Equations from strings.
    """

    # define abstract attributes from the Operator base.
    nin = 0
    nout = 1

    def __init__(self, name):
        Operator.__init__(self)
        ParameterSet.__init__(self, name)
        self.profile = None
        self.meta = {}
        return

    @property
    def symbol(self):
        """The symbolic representation, which is the generator name."""
        return self.name

    # Overload me!

    def __call__(self, x):
        """Evaluate the profile.

        This method must be overloaded. It only takes the independent
        variable to calculate over.

        Parameters
        ----------
        x : ndarray
            The independent variable over which to calculate.

        Returns
        -------
        ndarray
            The calculated profile.
        """
        return x

    # No need to overload anything below here

    def operation(self):
        """Evaluate the profile over the calculation points.

        Returns
        -------
        ndarray
            The profile calculated by ``__call__`` at ``profile.x``.
        """
        y = self.__call__(self.profile.x)
        return y

    def set_profile(self, profile):
        """Assign the profile.

        The generator stops observing any previously assigned Profile and
        starts observing the new one, so that changes to the Profile
        update ``meta`` and invalidate the cached value.

        Parameters
        ----------
        profile : Profile
            The Profile that specifies the calculation points and which
            will store the calculated signal.
        """
        if self.profile is not None:
            self.profile.remove_observer(self._on_profile_update)

        self.profile = profile
        self.profile.add_observer(self._on_profile_update)
        self._on_profile_update(other=(self,))
        return

    def _on_profile_update(self, other=()):
        if self.profile is not None:
            self.meta.update(self.profile.meta)
            self._process_metadata()
        self._flush(other=other)
        return

    def _process_metadata(self):
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

        Raises
        ------
        SrFitError
            If validation fails.
        """
        if self.profile is None:
            raise SrFitError("profile is None")
        self.profile._validate()
        ParameterSet._validate(self)
        return


# End class ProfileGenerator
