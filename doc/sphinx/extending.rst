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
into various places. For example, external forward calculators can be wrapped
up inside ProfileGenerators without modifying the calculator. Structure
adapters defined in the diffpy.srfit.structure module are also built around
this principle. These adapters are hierarchical ParameterSets (found in
diffpy.srfit.fitbase.parameterset) that encapsulate the different pieces of a
structure.  For example, the DiffpyStructure structure adapter in
diffpy.srfit.structure.diffpystructure contains a LatticeParSet that
encapsulates the lattice data and one AtomParSet per atom.  These each contain
parameters for what they encapsulate, such as lattice parameters or atom
positions. 

Fundamentally, it is the adjustable parameters of a structure container,
forward calculator or other object that needs to be adapted so that SrFit can
manipulate the underlying data object. These adapted parameters can then be
organized into ParameterSets, as in the case of a structure adapter. The
ParameterAdapter class found in diffpy.srfit.fitbase.parameter is designed for
this purpose. ParameterAdapter is a Parameter that defers to another object
when setting or retrieving its value.

We give a simplified example of how to do such an adaptation with a
hypothetical atom object called SimpleAtom that has attributes ``x``, ``y`` and
``z``. ::

    class SimpleAtomParSet(ParameterSet):
        
        def __init__(self, atom, name):
            ParameterSet.__init__(self, name)
            self.atom = atom

            # Here we create a Parameter for each x, y and z
            xpar = ParameterAdapter("x", atom, attr = "x")
            ypar = ParameterAdapter("y", atom, attr = "y")
            zpar = ParameterAdapter("z", atom, attr = "z")

            # Here we add them to the structure
            self.addParameter(xpar)
            self.addParameter(ypar)
            self.addParameter(zpar)

            return

ParameterAdapter as used here accesses the attribute ``x`` via ``atom.x`` by
specifing ``attr = "x"``. When retrieving or setting the value of ``xpar``, the
call would be redirected to ``atom.x``.

If SimpleAtom did not have an attribute named ``x``, but rather accessor
methods named ``getX`` and ``setX``, then the ParameterAdapter would be used
as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.getX, 
            setter = SimpleAtom.setX)

Note that the unbound methods are used. The names ``getter`` and ``setter``
describe how the accessor attributes are used to access the value of the
parameter. When ``xpar.getValue()`` is called, it redirects to
``SimpleAtom.getX(atom)``.

A third usage of ParameterAdapter is also possible. If instead SimpleAtom had a
methods called ``get`` and ``set`` that take as the second argument the name of
what to retrieve, then this can be adapted as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.get, 
            setter = SimpleAtom.set, attr = "x")

Thus, when ``xpar.getValue()`` is called, it in turn calls
``SimpleAtom.get(atom, "x")``. ``xpar.setValue(value)`` calls
``SimpleAtom.set(atom, "x", value)``.

If the attributes of an object object cannot be accessed in one of these three
ways, then you must write external accessor methods that can be set as the
getter and setter of the ParameterAdapter. For example, if the ``x``, ``y`` and
``z`` values were held in a list called ``xyz``, then you would have to write
the functions ``getX`` and ``setX`` that would manipulate this list, and use
these functions as in the second example.


Creating Custom FitHooks
--------------------------

The FitHook class is used by a FitRecipe to report fit progress to a user.
FitHook can be found in the diffpy.srfit.fitbase.fithook module.  FitHook can
be customized to provide customized fit output, such as a live plot of the
output. The FitHook class has three methods that one can overload.

.. currentmodule:: diffpy.srfit.fitbase.fithook

* .. automethod:: FitHook.reset
* .. automethod:: FitHook.precall
* .. automethod:: FitHook.postcall

To use a custom FitHook, assign an instance to a FitRecipe using the
``setFitHook`` method.


Extending Profile Parsers
--------------------------

The ProfileParser class is located in the diffpy.srfit.fitbase.parser module.
The purpose of this class is to read data and metadata from a file or string
and pass that data and metadata to a Profile instance. The Profile in turn will
pass this information to a ProfileGenerator.

The simplest way to extend the ProfileParser is to derive a new class from
ProfileParser and overload the parseString class. By default, the parseFile
class can read an ascii file and passes the loaded string to the parseString
method. For non-ascii data one should overload both of these methods. An
example of a customized ProfileParser is the PDFParser class in the
diffpy.srfit.pdf.pdfparser module.

Here is a simple example demonstrating how to extract (x,y) data from a
two-column string. ::

    def parseString(self, datastring):

        xvals = []
        yvals = []
        dxvals = None
        dyvals = None

        for line in datastring.split("\n"):
            sline = line.split()
            x = float(sline[0])
            y = float(sline[1])

            xvals.append(x)
            yvals.append(y)

        self._banks.append([xvals, yvals, dxvals, dyvals])
        return

The ``self._banks.append`` line puts the data arrays into the ``_banks`` list.
This list is for collecting multiple data sets that may be present within a
single file.  The ``dxvals`` and ``dyvals`` are the uncertainty values on the
``xvals`` and ``yvals``.  In this simple example they are not present, and so
are set to None.

In general, the data string may contain metadata. The ProfileParser has a
dictionary attribute named ``_meta``. The parser can put any information into
this dictionary. It is up to a ProfileGenerator that may use the parsed data
to look for usable information within this dictionary.

If the data is not in a form that can be used by a ProfileGenerator then
it is the responsibility of the parser to convert this data to a usable form.


Extending Profiles
--------------------------

Even with the ability to customize ProfileParsers, it may be necessary to
create custom Profile objects for different types of data. This is useful when
adapting an external data container to the SrFit interface. For example, the
external container may need to be retained so it can be used within an external
program before or after iterfacing with SrFit. An example of a customized
Profile is the SASProfile class in the diffpy.srfit.sas.sasprofile module:

.. literalinclude:: ../../diffpy/srfit/sas/sasprofile.py
   :pyobject: SASProfile

The ``__init__`` method sets the ``xobs``, ``yobs`` and ``dyobs`` attributes of
the SASProfile to the associated arrays within the ``_datainfo`` attribute. The
``setObservedProfile`` method is overloaded to modify the ``_datainfo`` arrays
when their corresponding attributes are modified. This keeps the arrays in
sync.


Custom Constraints and Restraints
----------------------------------------

**TODO**. For now, see the diffpy.srfit.structure.objcryststructure module.
