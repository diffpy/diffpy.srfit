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

#ifndef SRFIT_EVALUATOR_H
#define SRFIT_EVALUATOR_H

#include "Visitor.hpp"
#include "Clicker.hpp"
#include <boost/python.hpp>

namespace bp = boost::python;

namespace diffpy {
namespace srfit {

class Argument;
class Operator;

class Evaluator : public Visitor
{

    public:

    Evaluator();
    virtual ~Evaluator();

    virtual void visitArgument(Argument& a);
    virtual void visitOperator(Operator& o);

    Clicker clicker;

    bp::list vals;
    bp::object value;


}; 

} // namespace srfit
} // namespace diffpy

#endif // SRFIT_EVALUATOR_H

