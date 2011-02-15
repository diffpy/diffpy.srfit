#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Example of fitting a Gaussian to simulated data.

This is an example of building a fit recipe that can be driven by an optimizer
to fit a Gaussian profile to simulated data.  The purpose of this example is
familiarize developers and users with how to define a refinement with SrFit.

Instructions

If you have not yet run the example, run it and inspect the output. The example
script is driven by the 'main' method defined below. Take a look at that
function to learn how the refinement is built and executed.

Extensions

After reading through the code, try to perform the following tasks. The process
will leave you with a much better understanding of how SrFit works.

- Use scipy's 'fmin' optimizer to refine the residual. The 'res' object can be
  passed as a function for scalar optimizers.
- Modify the fitting equation by adding a constant. See if that improves the
  fit.
- Create a script that fits the Lorentzian profile in the data/lorentzian.dat
  file.

"""


# This imports everything a user needs to use SrFit.
from diffpy.srfit import *

def main():
    """Create the fit and optimize it to find the Gaussian parameters."""

    # To start, we need to load the data. The 'loadProfile' function can be
    # used for this. It returns a 'Profile' object that contains the data and
    # metadata, and can be used to resize the data arrays.
    profile = loadProfile("data/gaussian.dat")
    # The profile is iterable and returns parameters containing the data
    # arrays. 'x' is the independent variable, 'y' is the profile and 'dy' is
    # the uncertainty in the profile.
    x, y, dy = profile
    # These parameters act like references to the data arrays. They
    # have a 'name' attribute containing the name of the reference and a
    # 'value' attribute that holds the data values. We need these references to
    # construct the fitting equations below.
    print x.name
    print x.value
    # The data arrays can be extracted from the profile using the data
    # property.
    xdata, ydata, dydata = profile.data

    # Next we make variables for the fit. Variables are parameters, but they
    # have been flagged to be optimized.  There are a few ways to make
    # variables. We can create one variable using the 'Var' function.  The
    # function accepts the name and initial value for the variable.
    A = Var("A", 1.0)
    # We can also create multiple variables at once using 'Vars'. This function
    # accepts any number of name, value tuples - each one creates a variable.
    x0, sigma = Vars(("x0", 5), ("sigma", 1))

    # Next we symbolically create a fitting equation that simulates the data.
    # When everything was SrFit was imported, we also imported several standard
    # numpy functions that have been wrapped as symbolic operators. These
    # symbolic operators are named according to their numeric function with a
    # trailing underscore "_". We write the fitting equation with standard
    # python syntax using these functions and the variables and parameters
    # created above.
    fiteq = A * exp_(-0.5*(x-x0)**2/sigma**2)
    # The result is a symbolic function. It has the same interface as a
    # parameter. It has a name and it's value is computed from the parameters
    # and symbolic operations with the 'value' property.
    print fiteq.value

    # Now we have some data and the function that simulates the data so we can
    # create a residual object. The residual is what we will hand off to the
    # optimizer - the optimizer will try to minimize the residual function.  We
    # will do this the hard way to demonstrate the process. We want to optimize
    # the chi^2 function constructed from the fit equation and the data. We
    # need to construct the vector form of this, chi, where the scalar chi^2
    # can be computed from dot(chi, chi).
    chi = (fiteq - y) / dy
    # So the output looks nice, we rename this object.
    chi.name = "chi^2"
    # Now we create an object that can communicate more easily with an
    # optimizer.
    res = residual(chi)

    # Now we can refine the residual. The residual object we just created has a
    # 'vec' method that calculates the vector residual from 'chi'. (The object
    # can also be called as a function to compute the scalar residual.) The
    # 'leastsq' optimizer needs this function and a list of initial values to
    # work with.  We get these initial values from the 'values' property of
    # 'res'. 
    from scipy.optimize import leastsq
    leastsq(res.vec, res.values)

    # Now we have optimized the variables so that 'fiteq' best matches the
    # data. We can summarize the results with a FitResults object. This object
    # organizes the variables and computes their uncertainties and
    # correlations. We can 'show' and 'save' the results using this object.
    results = FitResults(res) 
    results.show()

    # We can also plot the results. We retrieve the numeric values of the
    # symbolic data and fit objects to do this.
    import pylab
    pylab.plot(x.value, y.value, 'bo')
    pylab.plot(x.value, fiteq.value, 'r-')
    pylab.show()

    return

if __name__ == "__main__":

    main()

# End of file
