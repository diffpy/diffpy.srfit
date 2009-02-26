#include "Literal.hpp"

using namespace diffpy::srfit;

Literal::Literal() {};

Literal::~Literal() {};

bp::object 
Literal::getValue()
{
    return this->value;
}
