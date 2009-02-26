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

#include "diffpy/srfit/UFuncOperator.hpp"
#include "diffpy/srfit/Operator.hpp"
#include "diffpy/srfit/Literal.hpp"

#include <boost/utility.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

namespace bp = boost::python;
using namespace boost::python;

using namespace diffpy::srfit;

namespace {

class UFuncOperatorWrap : public UFuncOperator,
                          public wrapper<UFuncOperator>
{

    public:


    UFuncOperatorWrap() : UFuncOperator() {};

    object default_callFunction(const object& o)
    {
        return UFuncOperator::callFunction(o);
    }

    object callFunction(const object& o)
    {
        if (override callFunction = this->get_override("callFunction")) 
            return callFunction(o);
        return default_callFunction(o);
    }

};

}

BOOST_PYTHON_MODULE(_ufuncoperator)
{
    import_array();
    import_ufunc();

    class_<UFuncOperatorWrap, boost::noncopyable, bases<Operator> >
        ("UFuncOperator")
        .def("setUFunc", &UFuncOperator::setUFunc,
            with_custodian_and_ward<1,2>())
        .def("callFunction", &UFuncOperator::callFunction, 
            &UFuncOperatorWrap::default_callFunction)
        ;

    // Derive the converter for from-python objects from this
    class_<UFuncOperator, bases<UFuncOperatorWrap> >("_UFuncOperator")
        ;
}
