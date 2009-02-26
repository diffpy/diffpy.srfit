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
* boost::python bindings to diffpy::srfit::Evaluator
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

#include "diffpy/srfit/Evaluator.hpp"
#include "diffpy/srfit/Visitor.hpp"

#include <boost/utility.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

namespace bp = boost::python;
using namespace boost::python;

using namespace diffpy::srfit;

namespace {

}

BOOST_PYTHON_MODULE(_evaluator)
{
    import_array();
    import_ufunc();

    class_<Evaluator, bases<Visitor> >("Evaluator")
        .def_readonly("clicker", &Evaluator::clicker)
        .def_readonly("value", &Evaluator::value)
        ;

}
