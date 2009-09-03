////////////////////////////////////////////////////////////////////////
//
// diffpy            by DANSE Diffraction group
//                   Simon J. L. Billinge
//                   (c) 2009 Trustees of the Columbia University
//                   in the City of New York.  All rights reserved.
//
// File coded by:    Chris Farrow
//
// See AUTHORS.txt for a list of people who contributed.
// See LICENSE.txt for license information.
//
// c++ implementation of Clicker class.
//
////////////////////////////////////////////////////////////////////////


//#define NDEBUG
#include <cassert>
#include <string>
#include <sstream>
#include <set>

#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>


using namespace boost::python;

namespace {

class Clicker
{

    public:

    Clicker() : state(0) {};

    /// Add a Clicker that observes this one.
    void addObserver(Clicker* other)
    {
        assert(other != NULL);
        observers.insert(other);
        update();
    }

    /// Add a Clicker to observe
    void addSubject(Clicker* other)
    {
        assert(other != NULL);
        other->addObserver(this);
    }

    /// Remove an observer
    ///
    /// This has no effect when the passed Clicker is not an observer
    void removeObserver(Clicker* other)
    {
        assert(other != NULL);
        it = observers.find(other);
        if(it != observers.end())
        {
            observers.erase( it );
        }
    }

    /// Remove a subject
    ///
    /// This has no effect when the passed Clicker is not an observer
    void removeSubject(Clicker* other)
    {
        assert(other != NULL);
        other->removeObserver(this);
    }

    /// Indicate if the passed Clicker is observing this Clicker
    bool hasObserver(Clicker* other)
    {
        assert(other != NULL);
        return observers.end() != observers.find(other);
    }

    /// Indicate if the passed Clicker is observed by this Clicker
    bool hasSubject(Clicker* other)
    {
        assert(other != NULL);
        return other->hasObserver(this);
    }

    /// Increment the Clicker global state, and update
    inline void click()
    {
        ++Clicker::numclicks;
        update();
    }

    /// Update the local state and that of the observers
    ///
    /// This sets the local state to the global state and updates all
    /// observers.
    inline void update()
    {
        state = Clicker::numclicks;
        for(it = observers.begin(); it != observers.end(); ++it)
        {
            (*it)->update();
        }
    }

    /// Compare the local state of two Clickers
    inline int __cmp__(Clicker* other)
    {
        assert(other != NULL);
        return this->state - other->state;
    }

    /// global state : local state
    std::string __str__()
    {
        std::stringstream s;
        s << state << ":" << Clicker::numclicks;
        return s.str();
    }


    /* data members */
    static unsigned long long numclicks;
    unsigned long long state;

    std::set<Clicker*> observers;
    std::set<Clicker*>::iterator it;


};

unsigned long long Clicker::numclicks = 0;

// Doc-strings

const char* clickerDoc =
"Clicker class for recording state changes.\n";

const char* addObserverDoc =
"Add a Clicker that observes this one.\n";

const char* addSubjectDoc =
"Add a Clicker to observe.\n";

const char* removeObserverDoc =
"Remove an observer.\n\n\
This has no effect when the passed Clicker is not an observer.\n";

const char* removeSubjectDoc =
"Remove an subject.\n\n\
This has no effect when the passed Clicker is not a subject.\n";

const char* hasObserverDoc =
"Indicate if the passed Clicker is observing this Clicker.\n";

const char* hasSubjectDoc =
"Indicate if the passed Clicker is observed by this Clicker.\n";

const char* clickDoc =
"Increment the Clicker global state, and update.\n";

const char* updateDoc =
"Update the local state and that of the observers.\n\n\
This sets the local state to the global state and updates all observers.\n";

const char* cmpDoc =
"Compare the local state of two Clickers.\n";

const char* strDoc =
"global state : local state.\n";

}

BOOST_PYTHON_MODULE(_clicker)
{


    class_<Clicker>("Clicker", clickerDoc)
        .def("addObserver", &Clicker::addObserver, 
            with_custodian_and_ward<1,2>(), addObserverDoc)
        .def("addSubject", &Clicker::addSubject, 
            with_custodian_and_ward<2,1>(), addSubjectDoc)
        .def("removeObserver", &Clicker::removeObserver, removeObserverDoc)
        .def("removeSubject", &Clicker::removeSubject, removeSubjectDoc)
        .def("hasObserver", &Clicker::hasObserver, hasObserverDoc)
        .def("hasSubject", &Clicker::hasSubject, hasSubjectDoc)
        .def("click", &Clicker::click, clickDoc)
        .def("update", &Clicker::update, updateDoc)
        .def("__cmp__", &Clicker::__cmp__, cmpDoc)
        .def("__str__", &Clicker::__str__, strDoc)
        ;

}
