// -*- C++ -*-
/*
* diffpy.srfit      by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2008 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
* 
* File coded by:    Chris Farrow
* 
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
* 
*/

#ifndef SRFIT_LITERAL_H
#define SRFIT_LITERAL_H

#include <boost/python.hpp>

#include <string>

#include "Clicker.hpp"

namespace bp = boost::python;

namespace diffpy {
namespace srfit {

class Visitor;

/* Abstract class for all visitors to a Literal tree. */
class Literal
{

    public:

    Literal();
    virtual ~Literal();

    virtual void identify(Visitor& v) = 0;

    bp::object getValue();

    // The name of this thing
    std::string name;
    // for holding the value
    bp::object value;
    // The clicker for the value
    Clicker clicker;

};
 
} // namespace srfit
} // namespace diffpy

#endif // SRFIT_LITERAL_H

