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

class ArgumentWrap : public Argument,
                    public wrapper<Argument>
{

    public:

    ArgumentWrap() : Argument() {};

    void default_identify(Visitor& v)
    {
        Argument::identify(v);
    }

    void identify(Visitor& v)
    {
        if (override identify = this->get_override("identify")) 
            identify(v);
        default_identify(v);
    }
};


}

BOOST_PYTHON_MODULE(_argument)
{
    import_array();
    import_ufunc();

    // Override the virtual functions in the public class
    class_<ArgumentWrap, boost::noncopyable ,bases<Literal> >("Argument")
        .def("identify", &Argument::identify, &ArgumentWrap::default_identify)
        .def("setValue", &Argument::setValue)
        ;

    // Derive the converter for from-python objects from this
    class_<Argument, bases<ArgumentWrap> >("_Argument");

}
