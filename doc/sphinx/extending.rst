.. _developers-guide-extending:

===================
Extending SrFit
===================

The examples in :ref:`_developers-guide-examples` give an overview of how to
use SrFit and extend it with various custom-made objects. Many pieces of SrFit
that are not covered in the examples are discussed here.

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
adapting an external data container to the SrFit interface. This is useful when
that external container can be used within an external program before or after
iterfacing with SrFit. An example of a customized Profile is the SASProfile
class in the diffpy.srfit.sas.sasprofile module.


Creating Custom FitHooks
--------------------------

The FitHook class is used by a FitRecipe to report fit progress to a user.
FitHook can be found in the diffpy.srfit.fitbase.fithook module.  FitHook can
be customized to provide customized fit output, such as a live plot of the
output. The FitHook class has three methods that one can overload.

 * ``reset`` - This is called whenever a FitRecipe is reconfigured after a
   change in the data structure.
 * ``precall`` - This gets called within FitRecipe.residual before the residual
   is calculated. It takes the FitRecipe instance as an argument.
 * ``postcall`` - This gets called within FitRecipe.residual after the residual
   is calculated. It takes the FitRecipe instance and residual vector as
   arguments.

To use a custom FitHook, assign an instance to the ``fithook`` attribute of a
FitRecipe.


Creating Custom Structure Containers
-------------------------------------

The structure containers for SrFit are found in the diffpy.srfit.structure
module. These containers are hierarchical ParameterSets (found in
diffpy.srfit.fitbase.parameterset) that encapsulate the sub-components of a
structure. For example, the DiffpyStructure structure adapter contains a
LatticeParSet that encapsulates the lattice data and one AtomParSet per atom.
These each contain parameters for the component, such as lattice parameters or
atom positions.

The purpose of a structure ParameterSet is to adapt a structure object so that
changing the Parameters in the ParameterSet mutates the structure object. The
mutated structure can then be used by a forward calculator that knows nothing
about SrFit. Just as well, changing the structure object directly (within
limitations) should reflect in the ParameterSet.

Useful in adapting an object in this way is the ParameterAdapter class found in
diffpy.srfit.fitbase.parameter. This class is designed to defer to another
object when setting or retrieving its value.

We give a simplified example of how to do such an adaptation with a
hypothetical atom object called SimpleAtom that has attributes x, y and z. ::

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

ParameterAdapter as used here accesses the attribute x via ``atom.x`` by
specifing ``attr = "x"``. If SimpleAtom did not have an attribute named x, but
rather accessor methods named ``getX`` and ``setX``, then the ParameterAdapter
would be used as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.getX, setter =
                SimpleAtom.setX)

Note that the unbound methods are used. A third usage of ParameterAdapter is
also possible. If instead SimpleAtom had a methods called ``get`` and ``set``
that take as the first argument the name of what to get, then this can be
adapted as::

    xpar = ParameterAdapter("x", atom, getter = SimpleAtom.get, setter =
                SimpleAtom.set, attr = "x")

If the attributes of the structure object cannot be accessed in one of these
three ways, then it is recommended to write external accessor methods that can
be set as the getter and setter of the ParameterAdapter. An example of how this
is done can be found in the LatticeParSet class in
diffpy.srfit.structure.diffpystructure.

Note that parameter adapters can be used wherever parameters can be. For
example, the built-in PDFGenerator in diffpy.srfit.pdf.pdfgenerator works in
concert with the structure adapters and an external PDF forward calculator by
adapting the parameters of the calculator. All of this is done without the
structure objects or forward calculator knowing anything about SrFit.


