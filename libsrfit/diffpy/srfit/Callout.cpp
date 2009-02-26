#include "Callout.hpp"
#include "Argument.hpp"
#include "Operator.hpp"
#include <vector>
#include <iostream>

using namespace diffpy::srfit;

Callout::Callout() {};

Callout::~Callout() {};

void 
Callout::visitArgument(Argument& a)
{
    std::cout << a.name << std::endl;
}

void 
Callout::visitOperator(Operator& o)
{
    std::cout << o.name << std::endl;
    for(std::vector<Literal*>::iterator lit = o.args.begin();
        lit != o.args.end(); ++lit)
    {

        (*lit)->identify(*this);
    }
}

