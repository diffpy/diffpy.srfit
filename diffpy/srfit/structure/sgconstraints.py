#!/usr/bin/env python
"""Code to set space group constraints for a crystal structure."""

__id__ = "$Id$"


def constrainSpaceGroup(lattice, system):
    """Constrain the lattice parameters according to the space group.

    This forces the lattice parameters to conform to the space group symmetry.
    The protocol this follows is listed under Crystal Systems below. 


    Arguments:

    lattice --  A ParameterSet representing a lattice. It must have the lattice
                parameters names "a", "b", "c", "alpha", "beta", "gamma".
    system  --  The crystal system (string), with one of the names below.
    
    Crystal Systems:

    Triclinic       --  No constraints.
    Monoclinic      --  alpha and beta are constrained to 90 unless alpha !=
                        beta and alpha == gamma, in which case alpha and
                        gamma are constrained to 90.
    Orthorhombic    --  alpha, beta and gamma are constrained to 90
    Tetragonal      --  b is constrained to a and alpha, beta and gamma are
                        constrained to 90.
    Trigonal        --  If gamma == 120, then b is constrained to a, alpha
                        and beta are constrained to 90 and gamma is
                        constrained to 120.  Otherwise, b and c are
                        constrained to a, beta and gamma are constrained to
                        alpha.
    Hexagonal       --  b is constrained to a, alpha and beta are
                        constrained to 90 and gamma is constrained to 120.
    Cubic           --  b and c are constrained to a, and alpha, beta and
                        gamma are constrained to 90.

    """
    # First clear any constraints or constant variables in the lattice
    lattice._constraints = {}
    for par in lattice._parameters:
        par.setConst(False)

    f = _constraintMap[system]
    f(lattice)
    return

def _constrainTriclinic(lattice):
    """Make constraints for Triclinic systems.

    This frees the current value of all parameters.
    """
    return

def _constrainMonoclinic(lattice):
    """Make constraints for Monoclinic systems.
    
    alpha and beta are constrained to 90 unless alpha != beta and alpha ==
    gamma, in which case alpha and gamma are constrained to 90.
    """
    lattice.alpha.setConst(True, 90.0)
    beta = lattice.beta.getValue()
    gamma = lattice.gamma.getValue()

    if 90 != beta and beta == gamma:
        lattice.gamma.setConst(True, 90)
    else:
        lattice.beta.setConst(True, 90)
    return

def _constrainOrthorhombic(lattice):
    """Make constraints for Orthorhombic systems.
    
    alpha, beta and gamma are constrained to 90
    """
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    return

def _constrainTetragonal(lattice):
    """Make constraints for Tetragonal systems.

    b is constrained to a and alpha, beta and gamma are constrained to 90.
    """
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    lattice.constrain(lattice.b, lattice.a)
    return

def _constrainTrigonal(lattice):
    """Make constraints for Trigonal systems.

    If gamma == 120, then b is constrained to a, alpha and beta are
    constrained to 90 and gamma is constrained to 120. Otherwise, b and c
    are constrained to a, beta and gamma are constrained to alpha.
    """
    if lattice.gamma.getValue() == 120:
        lattice.constrain(lattice.b, lattice.a)
        lattice.alpha.setConst(True, 90.0)
        lattice.beta.setConst(True, 90.0)
        lattice.gamma.setConst(True, 120)
    else:
        lattice.constrain(lattice.b, lattice.a)
        lattice.constrain(lattice.c, lattice.a)
        lattice.constrain(lattice.beta, lattice.alpha)
        lattice.constrain(lattice.gamma, lattice.alpha)
    return

def _constrainHexagonal(lattice):
    """Make constraints for Hexagonal systems.

    b is constrained to a, alpha and beta are constrained to 90 and gamma is
    constrained to 120.
    """
    lattice.constrain(lattice.b, lattice.a)
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 120.0)
    return

def _constrainCubic(lattice):
    """Make constraints for Cubic systems.

    b and c are constrained to a, alpha, beta and gamma are constrained to
    90.
    """
    lattice.constrain(lattice.b, lattice.a)
    lattice.constrain(lattice.c, lattice.a)
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    return

# This is used to map the correct crystal system to the proper constraint
# function. It is used by constrainSpaceGroup.
_constraintMap = {
  "Triclinic"  : _constrainTriclinic,
  "Monoclinic" : _constrainMonoclinic,
  "Orthorhombic" : _constrainOrthorhombic, 
  "Tetragonal" : _constrainTetragonal,
  "Trigonal"   : _constrainTrigonal,
  "Hexagonal"  : _constrainHexagonal,
  "Cubic"      : _constrainCubic
}


