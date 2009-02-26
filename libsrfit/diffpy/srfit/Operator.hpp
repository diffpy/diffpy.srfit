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

#ifndef SRFIT_OPERATOR_H
#define SRFIT_OPERATOR_H

#include "Literal.hpp"

#include <boost/python.hpp>
#include <string>
#include <vector>

namespace bp = boost::python;


namespace diffpy {
namespace srfit {

class Visitor;

class Operator : public Literal
{

    public:

    Operator();
    virtual ~Operator();

    virtual void identify(Visitor& v);
    void addLiteral(Literal& l);

    // This will let us keep a list of python types within this c++ object. This
    // is necessary to keep the identity of these objects for use in python.
    void addPyLiteral(bp::object& pyl);

    virtual bp::object callFunction(const bp::object&);

    std::string symbol;
    size_t nin;
    size_t nout;

    std::vector<Literal*> args;
    bp::list pyargs;

};
 
} // namespace srfit
} // namespace diffpy

#endif // SRFIT_OPERATOR_H

