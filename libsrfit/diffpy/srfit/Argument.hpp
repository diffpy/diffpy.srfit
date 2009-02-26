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

#ifndef SRFIT_ARGUMENT_H
#define SRFIT_ARGUMENT_H

#include "Literal.hpp"

#include <boost/python.hpp>

namespace bp = boost::python;

namespace diffpy {
namespace srfit {

class Argument : public Literal
{

    public:

    virtual ~Argument();

    void setValue(bp::object& val);

    virtual void identify(Visitor& v);

};
 
} // namespace srfit
} // namespace diffpy

#endif // SRFIT_ARGUMENT_H

