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
* boost::python bindings to diffpy::srfit::Operator
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

#include "diffpy/srfit/Operator.hpp"
#include "diffpy/srfit/Literal.hpp"
#include "diffpy/srfit/Visitor.hpp"

#include <boost/utility.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

namespace bp = boost::python;
using namespace boost::python;

using namespace diffpy::srfit;

namespace {

class OperatorWrap : public Operator,
                     public wrapper<Operator>
{

    public:


    OperatorWrap() : Operator() {};

    object default_callFunction(const object& o)
    {
        return Operator::callFunction(o);
    }

    object callFunction(const object& o)
    {
        if (override callFunction = this->get_override("callFunction")) 
            return callFunction(o);
        return default_callFunction(o);
    }

    void default_identify(Visitor& v)
    {
        Operator::identify(v);
    }

    void identify(Visitor& v)
    {
        if (override identify = this->get_override("identify")) 
            identify(v);
        default_identify(v);
    }

};

}

BOOST_PYTHON_MODULE(_operator)
{
    import_array();
    import_ufunc();
    // Override the virtual functions in the public class
    class_<OperatorWrap, boost::noncopyable, bases<Literal> >("Operator")
        .def("identify", &Operator::identify, &OperatorWrap::default_identify)
        .def("callFunction", &Operator::callFunction, 
            &OperatorWrap::default_callFunction)
        .def("addLiteral", &OperatorWrap::addPyLiteral,
            with_custodian_and_ward<1,2>())
        .def_readwrite("nin", &Operator::nin)
        .def_readwrite("nout", &Operator::nout)
        .def_readwrite("symbol", &Operator::symbol)
        .def_readonly("args", &Operator::pyargs)
        ;

    // Derive the converter for from-python objects from this
    class_<Operator, bases<OperatorWrap> >("_Operator")
        ;


}
