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
* boost::python bindings to diffpy::srfit::Literal
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

#include "diffpy/srfit/Argument.hpp"
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

size_t getNumClicks()
{
    return Clicker::numclicks;
}

}

BOOST_PYTHON_MODULE(_clicker)
{
    import_array();
    import_ufunc();

    // Override the virtual functions in the public class
    class_<Clicker>("Clicker")
        .def("addObserver", &Clicker::addObserver,
            with_custodian_and_ward<1,2>())
        .def("addSubject", &Clicker::addSubject,
            with_custodian_and_ward<2,1>())
        .def("removeObserver", &Clicker::removeObserver)
        .def("removeSubject", &Clicker::removeSubject)
        .def("hasObserver", &Clicker::hasObserver)
        .def("hasSubject", &Clicker::hasSubject)
        .def("click", &Clicker::click)
        .def("update", &Clicker::update)
        .def(self < self)
        .def(self <= self)
        .def(self > self)
        .def(self >= self)
        .def(self == self)
        .def_readonly("state", &Clicker::state)
        .add_static_property("numclicks", &getNumClicks)
        ;



}
