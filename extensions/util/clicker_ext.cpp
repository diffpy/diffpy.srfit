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

typedef unsigned long long state_type;

class Clicker
{

    public:

    Clicker() : state(0) {};

    /// Add a Clicker that observes this one.
    void addObserver(Clicker* other)
    {
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

        observers.insert(other);
        update();
    }

    /// Add a Clicker to observe
    void addSubject(Clicker* other)
    {
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

        other->addObserver(this);
    }

    /// Remove an observer
    ///
    /// This has no effect when the passed Clicker is not an observer
    void removeObserver(Clicker* other)
    {
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

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
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

        other->removeObserver(this);
    }

    /// Indicate if the passed Clicker is observing this Clicker
    bool hasObserver(Clicker* other)
    {
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

        return observers.end() != observers.find(other);
    }

    /// Indicate if the passed Clicker is observed by this Clicker
    bool hasSubject(Clicker* other)
    {
        if( other == NULL )
        {
            PyErr_SetString(PyExc_ValueError, "Cannot add null clicker");
            throw_error_already_set();
        }

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
    inline int __cmp__(Clicker& other)
    {
        return state - other.state;
    }

    /// global state : local state
    std::string __str__()
    {
        std::stringstream s;
        s << state << ":" << Clicker::numclicks;
        return s.str();
    }


    /* data members */
    static state_type numclicks;
    state_type state;

    std::set<Clicker*> observers;
    std::set<Clicker*>::iterator it;

};

state_type Clicker::numclicks = 0;


// FlatClicker - This clicker class avoids walking the network during update.
// FlatClicker has extra overhead in informing all observers when a new clicker
// is added to or removed from the network.

class FlatClicker
{

    public:

    FlatClicker() : state(0) {};

    /// Add the clicker as an observer all the way down the tree.
    /// Do not expose
    void addObserverDown(FlatClicker* other)
    {
        observers.insert(other);

        // Pass the observer on to subjects
        std::set<FlatClicker*>::iterator it;
        for(it=subjects.begin(); it!=subjects.end(); ++it)
        {
            (*it)->addObserverDown(other);
        }
        return;
    }

    /// Add a FlatClicker that observes this one.
    void addObserver(FlatClicker& other)
    {
        // Make sure we don't create loops
        if( &other == this ||
            other.observers.find(this) != other.observers.end())
        {
            PyErr_SetString(PyExc_ValueError, "Cannot create loop");
            throw_error_already_set();
        }

        // Add the new observer and pass it down.
        addObserverDown(&other);
        // Make we are a direct subject of the observer.
        other.subjects.insert(this);
        // Add other's observers as well
        std::set<FlatClicker*>::iterator it;
        for(it = other.observers.begin(); it != other.observers.end(); ++it)
        {
            addObserverDown(*it);
        }
        update();
    }

    /// Add a FlatClicker to observe
    void addSubject(FlatClicker& other)
    {
        other.addObserver(*this);
    }

    /// Remove the clicker as an observer all the way down the tree.
    /// Does not remove the clicker at the top level.
    /// Do not expose
    void removeObserverDown(FlatClicker& other)
    {
        // If we're one of other's direct subject's, then we don't want to
        // break the link. We only break the link in removeObserver.
        if( other.subjects.find(this) != other.subjects.end() )
        {
            return;
        }

        observers.erase( &other );

        std::set<FlatClicker*>::iterator it;
        for(it=subjects.begin(); it!=subjects.end(); ++it)
        {
            (*it)->removeObserverDown(other);
        }
    }

    /// Remove an observer
    ///
    /// This has no effect when the passed FlatClicker is not an observer
    void removeObserver(FlatClicker& other)
    {
        if(observers.find(&other) != observers.end())
        {
            other.subjects.erase(this);
            removeObserverDown(other);
        }
    }

    /// Remove a subject
    ///
    /// This has no effect when the passed FlatClicker is not an observer
    void removeSubject(FlatClicker& other)
    {
        other.removeObserver(*this);
    }

    /// Indicate if the passed FlatClicker is a direct observer of this FlatClicker
    bool hasObserver(FlatClicker& other)
    {
        return other.hasSubject(*this);
    }

    /// Indicate if the passed FlatClicker is observed by this FlatClicker
    bool hasSubject(FlatClicker& other)
    {
        return subjects.end() != subjects.find(&other);
    }

    /// Increment the FlatClicker global state, and update
    inline void click()
    {
        ++FlatClicker::numclicks;
        update();
    }

    /// Update the local state and that of the observers
    ///
    /// This sets the local state to the global state and updates all
    /// observers.
    inline void update()
    {
        state = FlatClicker::numclicks;
        for(it = observers.begin(); it != observers.end(); ++it)
        {
            (*it)->state = state;
        }
    }

    /// Compare the local state of two FlatClickers
    inline int __cmp__(FlatClicker& other)
    {
        return state - other.state;
    }

    /// global state : local state
    std::string __str__()
    {
        std::stringstream s;
        s << state << ":" << FlatClicker::numclicks;
        return s.str();
    }


    /* data members */
    static state_type numclicks;
    state_type state;

    // All observers, throughout the whole network
    std::set<FlatClicker*> observers;
    // Just our direct subjects
    std::set<FlatClicker*> subjects;
    std::set<FlatClicker*>::iterator it;

};

state_type FlatClicker::numclicks = 0;

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

// Factor this out so we can wrap multiple clicker classes.
template <class C>
void
wrapClickerClass(class_<C>& c)
{
        c.def("addObserver", &C::addObserver, 
            with_custodian_and_ward<1,2>(), addObserverDoc)
        .def("addSubject", &C::addSubject, 
            with_custodian_and_ward<2,1>(), addSubjectDoc)
        .def("removeObserver", &C::removeObserver, removeObserverDoc)
        .def("removeSubject", &C::removeSubject, removeSubjectDoc)
        .def("hasObserver", &C::hasObserver, hasObserverDoc)
        .def("hasSubject", &C::hasSubject, hasSubjectDoc)
        .def("click", &C::click, clickDoc)
        .def("update", &C::update, updateDoc)
        .def("__cmp__", &C::__cmp__, cmpDoc)
        .def("__str__", &C::__str__, strDoc)
        ;
}


BOOST_PYTHON_MODULE(_clicker)
{

    class_<Clicker> PythonClicker("Clicker", clickerDoc);
    wrapClickerClass<Clicker>(PythonClicker);

    class_<FlatClicker> PythonFlatClicker("FlatClicker", clickerDoc);
    wrapClickerClass<FlatClicker>(PythonFlatClicker);
}
