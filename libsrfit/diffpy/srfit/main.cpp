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

#define PY_ARRAY_UNIQUE_SYMBOL _SRFIT_PY_ARRAY
#define PY_UFUNC_UNIQUE_SYMBOL _SRFIT_PY_UFUNC

#include <boost/python.hpp>

#include <numpy/arrayobject.h>
#include <numpy/noprefix.h>
#include <numpy/ufuncobject.h>

#include <string>
#include <iostream>

#include "UFuncOperator.hpp"
#include "Argument.hpp"
#include "Callout.hpp"

using namespace std;
using namespace boost::python;

using namespace diffpy::srfit;

void testCallout()
{

    UFuncOperator o;
    o.name = std::string("operator");
    Argument arg1;
    arg1.name = std::string("arg1");
    Argument arg2;
    arg2.name = std::string("arg2");

    o.addLiteral(arg1);
    o.addLiteral(arg2);

    Callout c;

    o.identify(c);
}


void testDefaultOperator()
{
    Operator o;
    object args(1.1);

    object retval = o.callFunction(args);
    cout << extract<float>(retval) << endl;

}

int main()
{
    testCallout();
    testDefaultOperator();

    return 0;
}

