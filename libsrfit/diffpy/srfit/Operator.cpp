#include "Operator.hpp"
#include "Literal.hpp"
#include "Visitor.hpp"

#include <boost/python.hpp>

#include <string>
#include <iostream>

using namespace diffpy::srfit;
namespace bp = boost::python;

Operator::Operator() : Literal(), nin(2), nout(1) {};

Operator::~Operator() {};

void 
Operator::identify(Visitor& v)
{
    v.visitOperator(*this);
}

bp::object
Operator::callFunction(const bp::object& o)
{
    return o;
}

void
Operator::addPyLiteral(bp::object& pyl)
{
    pyargs.append(pyl);
    Literal& l = bp::extract<Literal&>(pyl);
    addLiteral(l);
}

void 
Operator::addLiteral(Literal& l)
{
    args.push_back(&l);
    clicker.addSubject(&(l.clicker));
}

