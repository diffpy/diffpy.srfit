#include "Argument.hpp"
#include "Visitor.hpp"

#include <boost/python.hpp>
#include <numpy/arrayobject.h>

using namespace diffpy::srfit;
namespace bp = boost::python;

Argument::~Argument() {};

void
Argument::setValue(bp::object& val)
{
    this->value = val;
    this->clicker.click();
}

void
Argument::identify(Visitor& v)
{
    v.visitArgument(*this);
}

