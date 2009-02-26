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

#ifndef SRFIT_CLICKER_H
#define SRFIT_CLICKER_H

#include <set>

namespace diffpy {
namespace srfit {

class Clicker
{

    public:

    Clicker();
    virtual ~Clicker();

    void addObserver(Clicker* c);
    void addSubject(Clicker* c);
    void removeObserver(Clicker* c);
    void removeSubject(Clicker* c);
    bool hasObserver(const Clicker* c) const;
    bool hasSubject(const Clicker* c) const;
    void click();
    void update();

    bool operator<(const Clicker& c) const;
    bool operator<=(const Clicker& c) const;
    bool operator>(const Clicker& c) const;
    bool operator>=(const Clicker& c) const;
    bool operator==(const Clicker& c) const;

    static size_t numclicks;
    size_t state;
    private:
    std::set<Clicker*> observers;

}; 


} // namespace srfit
} // namespace diffpy

#endif // SRFIT_CLICKER_H

