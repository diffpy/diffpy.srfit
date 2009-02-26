/*****************************************************************************
*
* pyobjcryst        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* boost::python bindings to diffpy::srfit::Visitor
*
*
* $Id$
*
*****************************************************************************/
#define PY_ARRAY_UNIQUE_SYMBOL _SRFIT_PY_ARRAY
#define PY_UFUNC_UNIQUE_SYMBOL _SRFIT_PY_UFUNC

#include <boost/python.hpp>

#include <numpy/arrayobject.h>
#include <numpy/noprefix.h>
#include <numpy/ufuncobject.h>

#include "diffpy/srfit/Visitor.hpp"
#include "diffpy/srfit/Argument.hpp"
#include "diffpy/srfit/Operator.hpp"

#include <boost/utility.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <iostream>

namespace bp = boost::python;
using namespace boost::python;

using namespace diffpy::srfit;

namespace {

class VisitorWrap : public Visitor,
                    public wrapper<Visitor>
{

    public:

    VisitorWrap() : Visitor() {};

    void default_visitArgument(Argument& l)
    {
        Visitor::visitArgument(l);
    }

    void visitArgument(Argument& l)
    {
        if (override visitArgument = this->get_override("visitArgument")) 
            visitArgument(l);
        default_visitArgument(l);
    }

    void default_visitOperator(Operator& l)
    {
        Visitor::visitOperator(l);
    }

    void visitOperator(Operator& l)
    {
        if (override visitOperator = this->get_override("visitOperator")) 
            visitOperator(l);
        default_visitOperator(l);
    }

};

}

BOOST_PYTHON_MODULE(_visitor)
{
    import_array();
    import_ufunc();

    // Override the virtual functions in the public class
    class_<VisitorWrap, boost::noncopyable >("Visitor")
        .def("visitArgument", &Visitor::visitArgument,
            &VisitorWrap::default_visitArgument)
        .def("visitOperator", &Visitor::visitOperator,
            &VisitorWrap::default_visitOperator)
        ;
    // Derive the converter for from-python objects from this
    class_<Visitor, bases<VisitorWrap> >("_Visitor")
        ;
}
