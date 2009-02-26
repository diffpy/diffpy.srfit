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
#include "Clicker.hpp"
#include <set>

using namespace diffpy::srfit;

size_t Clicker::numclicks = 0;

Clicker::Clicker() : state(0) {};

Clicker::~Clicker() {};

void Clicker::addObserver(Clicker* c)
{
    observers.insert(c);
    update();
}

void Clicker::addSubject(Clicker* c)
{
    c->addObserver(this);
}

void Clicker::removeObserver(Clicker* c)
{
    std::set<Clicker*>::iterator item = observers.find(c);
    if( item != observers.end() )
    {
        observers.erase(item);
    }
}

void Clicker::removeSubject(Clicker* c)
{
    c->removeObserver(this);
}

bool Clicker::hasObserver(const Clicker* c) const
{
    std::set<Clicker*>::const_iterator item; 
    item = observers.find(const_cast<Clicker*>(c));
    return (item != observers.end());
}

bool Clicker::hasSubject(const Clicker* c) const
{
    return c->hasObserver(this);
}

void Clicker::click()
{
    ++Clicker::numclicks;
    update();
}

void Clicker::update()
{
    state = Clicker::numclicks;
    std::set<Clicker*>::iterator item;
    for(item = observers.begin(); item != observers.end(); ++item)
    {
        (*item)->update();
    }
}

bool Clicker::operator<(const Clicker &c) const
{
    return (state < c.state);
}

bool Clicker::operator<=(const Clicker &c) const
{
    return (state <= c.state);
}

bool Clicker::operator>(const Clicker &c) const
{
    return (state > c.state);
}

bool Clicker::operator>=(const Clicker &c) const
{
    return (state >= c.state);
}

bool Clicker::operator==(const Clicker &c) const
{
    return (state == c.state);
}

