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

#ifndef SRFIT_VISITOR_H
#define SRFIT_VISITOR_H

namespace diffpy {
namespace srfit {

class Argument;
class Operator;

class Visitor
{

    public:

    Visitor();
    virtual ~Visitor();

    virtual void visitArgument(Argument& a);
    virtual void visitOperator(Operator& o);


}; 

} // namespace srfit
} // namespace diffpy

#endif // SRFIT_VISITOR_H

