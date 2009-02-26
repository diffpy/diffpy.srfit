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

class LiteralWrap : public Literal,
                    public wrapper<Literal>
{

    public:

    LiteralWrap() : Literal() {};

    void identify(Visitor& v)
    {
        this->get_override("identify")(v);
    }
};

}

BOOST_PYTHON_MODULE(_literal)
{

    import_array();
    import_ufunc();

    class_<LiteralWrap, boost::noncopyable >("Literal", no_init)
        .def("identify", pure_virtual(&LiteralWrap::identify))
        .def("getValue", &Literal::getValue)
        .def_readwrite("name", &Literal::name)
        .def_readonly("clicker", &Literal::clicker)
        ;

}
