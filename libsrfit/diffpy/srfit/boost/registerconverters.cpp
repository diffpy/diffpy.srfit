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
* boost::python bindings to ObjCryst::AsymmetricUnit.
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

#include <boost/python/module.hpp>


BOOST_PYTHON_MODULE(_registerconverters)
{
    // There is a problem with getting this to work with libsrfit.so
    import_array();
    import_ufunc();

}
