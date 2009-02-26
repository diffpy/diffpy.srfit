#include "Visitor.hpp"
#include "Argument.hpp"
#include "Operator.hpp"
#include <vector>
#include <iostream>

using namespace diffpy::srfit;

Visitor::Visitor() {};

Visitor::~Visitor() {};

void Visitor::visitArgument(Argument& a) {}

void Visitor::visitOperator(Operator& o) {}
