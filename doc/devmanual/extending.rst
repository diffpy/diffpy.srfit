.. _developers-guide-extending:

===================
Extending SrFit
===================

The :ref:`developers-guide-examples` give an overview of how to
use SrFit and extend it with various custom-made objects. Many pieces of SrFit
that are not covered in the examples are discussed here.

Plugging Other Objects into SrFit
-------------------------------------

Much of the power of SrFit comes from being able to plug existing python codes
into the framework. For example, external forward calculators can be wrapped up
inside ``ProfileGenerators`` without modifying the calculator. This is
demonstrated in :ref:`developers-guide-examples`. Structure adapters defined in
the diffpy.srfit.structure module are also built around this principle. These
adapters are hierarchical ``ParameterSets`` (found in
``diffpy.srfit.fitbase.parameterset``) that encapsulate the different pieces of
a structure.  For example, the ``DiffpyStructureParSet`` structure adapter in
``diffpy.srfit.structure.diffpyparset`` contains ``DiffpyLatticeParSet``, which
encapsulates the lattice data and one ``DiffpyAtomParSet`` per atom.  These
each contain parameters for what they encapsulate, such as lattice parameters
or atom positions. 

Fundamentally, it is the adjustable parameters of a structure container,
forward calculator or other object that needs to be adapted so that SrFit can
manipulate the underlying data object. These adapted parameters can then be
organized into ``ParameterSets``, as in the case of a structure adapter. The
``ParameterAdapter`` class found in ``diffpy.srfit.fitbase.parameter`` is
designed for this purpose. ``ParameterAdapter`` is a ``Parameter`` that defers
to another object when setting or retrieving its value.

.. currentmodule:: diffpy.srfit.fitbase.parameter
.. autoclass:: ParameterAdapter

    The `name` argument is used to give attribute access to the
    ParameterAdapter instance when it is added to a ParameterSet or similar
    object. The `obj` argument is the parameter-like object to be adapted. It
    must provide some form of access to its data. If it provides a getter and
    setter, these can be specified with the `getter` and `setter` arguments.
    If the getter and setter require an attribute name, this is specified with
    the `attr` argument. If the data can be retrieved as an attribute, then the
    name of this attribute can be passed in the `attr` argument.


Here is a simple example of using ``ParameterAdapter`` to adapt a hypothetical
atom object called ``SimpleAtom`` that has attributes ``x``, ``y`` and ``z``. 
::

    class SimpleAtom(object):
        """Simple class holding x, y and z coordinates of an atom."""

        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z
            return

    # End class SimpleAtom

    class SimpleAtomParSet(ParameterSet):
        """Class adapting the x, y and z attributes of SimpleAtom as Parameters.""" 

        def __init__(self, atom, name):
            ParameterSet.__init__(self, name)
            # Store the atom, we might need it later
            self.atom = atom

            # Create a ParameterAdapter for the x, y and z attributes of atom
            xpar = ParameterAdapter("x", atom, attr = "x")
            ypar = ParameterAdapter("y", atom, attr = "y")
            zpar = ParameterAdapter("z", atom, attr = "z")

            # Add these to the parameter set
            self.addParameter(xpar)
            self.addParameter(ypar)
            self.addParameter(zpar)

            return

    # End class SimpleAtomParSet

The ``x``, ``y`` and ``z`` attributes (specified by the ``attr`` keyword
argument of ``ParameterAdapter``) of a ``SimpleAtom`` are wrapped as
``ParameterAdapter`` objects named `x`, `y`, and `z`.  They are then added to
the ``SimpleAtomParSet`` using the ``addParameter`` method, which makes them
accessible as attributes. 

If SimpleAtom did not have an attribute named ``x``, but rather accessor
methods named ``getX`` and ``setX``, then the ``ParameterAdapter`` would be
used as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.getX, 
        setter = SimpleAtom.setX)

Note that the *unbound* methods are used. The names ``getter`` and ``setter``
describe how the accessor attributes are used to access the value of the
parameter. When ``xpar.getValue()`` is called, it redirects to
``SimpleAtom.getX(atom)``. 

If instead ``SimpleAtom`` had methods called ``get`` and ``set`` that take as
the second argument the name of the attribute to retrieve or modify, then this
can be adapted as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.get, 
            setter = SimpleAtom.set, attr = "x")

Thus, when ``xpar.getValue()`` is called, it in turn calls
``SimpleAtom.get(atom, "x")``. ``xpar.setValue(value)`` calls
``SimpleAtom.set(atom, "x", value)``.

If the attributes of an object cannot be accessed in one of these three ways,
then you must write external accessor methods that can be set as the getter and
setter of the ``ParameterAdapter``. For example, if the ``x``, ``y`` and ``z``
values were held in a list called ``xyz``, then you would have to write the
functions ``getX`` and ``setX`` that would manipulate this list, and use these
functions as in the second example.


Extending Profile Parsers
--------------------------

The ``ProfileParser`` class is located in the ``diffpy.srfit.fitbase.parser``
module.  The purpose of this class is to read data and metadata from a file or
string and pass those data and metadata to a ``Profile`` instance. The
``Profile`` in turn will pass this information to a ``ProfileGenerator``.

The simplest way to extend the ``ProfileParser`` is to derive a new class from
``ProfileParser`` and overload the ``parseString`` method. By default, the
``parseFile`` method can read an ASCII file and passes the loaded string to the
``parseString`` method. For non-ASCII data one should overload both of these
methods. An example of a customized ``ProfileParser`` is the ``PDFParser``
class in the ``diffpy.srfit.pdf.pdfparser`` module.

Here is a simple example demonstrating how to extract (x,y) data from a
two-column string. ::

    def parseString(self, datastring):

        xvals = []
        yvals = []
        dxvals = None
        dyvals = None

        for line in datastring.splitlines():

            sline = line.split()
            x, y = map(float, sline)
            xvals.append(x)
            yvals.append(y)

        self._banks.append([xvals, yvals, dxvals, dyvals])
        return

The ``self._banks.append`` line puts the data arrays into the ``_banks`` list.
This list is for collecting multiple data sets that may be present within a
single file.  The ``dxvals`` and ``dyvals`` are the uncertainty values on the
``xvals`` and ``yvals``.  In this simple example they are not present, and so
are set to None.

In general, the data string may contain metadata. The ``ProfileParser`` has a
dictionary attribute named ``_meta``. The parser can put any information into
this dictionary. It is up to a ``ProfileGenerator`` that may use the parsed
data to define and retrieve usable metadata.

If the data are not in a form that can be stored in a ``Profile`` then it is
the responsibility of the parser to convert this data to a usable form.


Extending Profiles
--------------------------

Even with the ability to customize ``ProfileParsers,`` it may be necessary to
create custom ``Profile`` objects for different types of data. This is useful
when adapting an external data container to the SrFit interface. For example,
the external container may need to be retained so it can be used within an
external program before or after interfacing with SrFit. An example of a
customized Profile is the ``SASProfile`` class in the
``diffpy.srfit.sas.sasprofile`` module:

.. literalinclude:: ../../diffpy/srfit/sas/sasprofile.py
   :pyobject: SASProfile

The ``__init__`` method sets the ``xobs``, ``yobs`` and ``dyobs`` attributes of
the ``SASProfile`` to the associated arrays within the ``_datainfo`` attribute.
The ``setObservedProfile`` method is overloaded to modify the ``_datainfo``
arrays when their corresponding attributes are modified. This keeps the arrays
in sync.


Custom Restraints 
----------------------------------------

Restraints in SrFit are one way to include known information about a system
into a fit recipe. When customizing SrFit for a specific purpose, one may want
to create restraints. One example of this is in the ``SrRealParSet`` base class
in ``diffpy.srfit.structure.srrealparset``. SrReal provides many real-space
structure utilities for compatible structures, such as a PDF calculator and a
bond-valence sum (BVS) calculator. The PDF calculator works very well as a
``ProfileGenerator`` (see the :ref:`developers-guide-examples`), but the BVS
calculator is better suited as a restraint. This makes it very easy to keep the
BVS constrained during a PDF fit or some other refinement.

Creating a custom restraint is a two-step process. First, a class must be
derived from ``diffpy.srfit.fitbase.restraint.Restraint`` that can calculate
the restraint cost.  This requires the ``penalty`` method to be overloaded.
This method has the following signature

.. currentmodule:: diffpy.srfit.fitbase.restraint
.. automethod:: Restraint.penalty

The `w` factor is optionally used to scale the restraint cost. Its purpose is
to keep the restraint cost comparable to the residual of a single data point.

``BVSRestraint`` from ``diffpy.srfit.structure.bvsrestraint`` is a custom
``Restraint`` whose penalty is the root-mean-square deviation from the expected
and calculated BVS of a structure.

.. literalinclude:: ../../diffpy/srfit/structure/bvsrestraint.py
   :pyobject: BVSRestraint

Note that the penalty scaling is optional (selected by the `scaled` flag) and
uncertainty on the result (`sig`) may be applied. These two options are
recommended with any custom ``Restraint``. 

The second part of a custom restraint is to allow it to be created from a
restrainable object. A ``BVSRestraint`` is used to restrain a ``SrRealParSet``,
which is a ``ParameterSet`` wrapper base class for SrReal-compatible
structures.  The restraint is applied with the ``restrainBVS`` method.

.. literalinclude:: ../../diffpy/srfit/structure/srrealparset.py
   :pyobject: SrRealParSet.restrainBVS

The purpose of the method is to create the custom ``Restraint`` object,
configure it and store it. Note that the optional `sig` and `scaled` flag are
passed as part of this method. Both ``_restraints`` and
``_updateConfiguration`` come from ``ParameterSet``, from which
``SrRealParSet`` is derived. The ``_restraints`` attribute is a set of
``Restraints`` on the object. The ``_updateConfiguration`` method makes any
object containing the ``SrRealParSet`` aware of the configuration change.  This
gets propagated to the top-level ``FitRecipe``, if there is one.  The restraint
object is returned by the method so that it may be later removed.

For more examples of custom restraints can be found in the
``diffpy.srfit.structure.objcrystparset`` module.


Custom FitHooks
--------------------------

The ``FitHook`` class is used by a ``FitRecipe`` to report fit progress to a
user.  ``FitHook`` can be found in the ``diffpy.srfit.fitbase.fithook`` module.
``FitHook`` can be customized to provide customized fit output, such as a live
plot of the output. The ``FitHook`` class has three methods that one can
overload.

.. currentmodule:: diffpy.srfit.fitbase.fithook

* .. automethod:: FitHook.reset
* .. automethod:: FitHook.precall
* .. automethod:: FitHook.postcall

To use a custom ``FitHook``, assign an instance to a ``FitRecipe`` using the
``pushFitHook`` method. All ``FitHook`` instances held by a ``FitRecipe`` will
be used in sequence during a call to ``FitRecipe.residual``.




