#include "Evaluator.hpp"
#include "Literal.hpp"
#include "Argument.hpp"
#include "Operator.hpp"

#include <vector>
#include <iostream>

#include <boost/python.hpp>

namespace bp = boost::python;

using namespace diffpy::srfit;

Evaluator::Evaluator() {};

Evaluator::~Evaluator() {};

void 
Evaluator::visitArgument(Argument& a)
{
    this->value = a.value;
    this->vals = bp::list();
    this->vals.append(a.value);
}

void 
Evaluator::visitOperator(Operator& o)
{
    if( o.clicker > this->clicker )
    {

        bp::list localvals;
        for(std::vector<Literal*>::iterator lit = o.args.begin(); 
            lit != o.args.end(); ++lit)
        {
            (*lit)->identify(*this);
            localvals.extend(this->vals);
        }

        //std::cout << "localvals " << bp::len(localvals) << std::endl;
        o.value = o.callFunction(bp::tuple(localvals));
        //std::cout << bp::extract<float>(o.value) << std::endl;
    }

    this->value = o.value;

    if( o.nout == 1 )
    {
        this->vals = bp::list();
        this->vals.append(this->value);
    }
    else
    {
        this->vals = bp::list(this->value);
    }
}

