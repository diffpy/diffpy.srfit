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

#ifndef SRFIT_UFUNCOPERATOR_H
#define SRFIT_UFUNCOPERATOR_H

#include "Operator.hpp"

namespace bp = boost::python;

namespace diffpy {
namespace srfit {

class Visitor;

class UFuncOperator : public Operator
{

    public:

    UFuncOperator();
    virtual ~UFuncOperator();

    virtual bp::object callFunction(const bp::object&);

    void setUFunc(bp::object&, std::string symbol = "");

    private:
    PyObject* f;

};

} // namespace srfit
} // namespace diffpy

#endif // SRFIT_UFUNCOPERATOR_H

