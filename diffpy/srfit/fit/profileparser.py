#!/usr/bin/env python
"""This module contains classes for parsing profiles from files.

ProfileParser is a base class for parsing data. It can interact with a Profile
object to automatically set the Profile's data and metadata. Each specific file
format must be encapsulated in a ProfileParser subclass.

See the class documentation for more information.

"""
__all__ = ["getParser", "parserInfo"]

# Registry of parsers indexed by parser format.
_registry = {}

def getParser(fmt):
    """Get a parser class based on its format.

    Raises ValueError if a parser for the format cannot be found.

    """
    try:
        ParserClass = _registry[fmt]
    except KeyError:
        msg = "Parser for '%s' format cannot be found"%fmt
        raise ValueError(msg)

    return ParserClass

def parserInfo():
    """Print information on all parsers."""
    import textwrap
    def _s(docstr):
        lines = docstr.splitlines()
        for line in lines:
            line.strip()
            if line: return line
        return ""
    maxlen = max(map(len, _registry.keys()))
    pad = maxlen + 2
    indent = pad + 4
    items = ((fmt, _s(cls.__doc__)) for fmt, cls in _registry.iteritems())
    lines = (fmt.ljust(pad) + "--  " + doc for fmt, doc in items)
    tw = textwrap.TextWrapper(width = 79, subsequent_indent = ' '*indent)
    wlines = (tw.fill(line) for line in lines)
    info = "Registered parsers\n\n"
    info += "\n".join(sorted(wlines))
    print info
    return

class ProfileParser(object):
    """Class for parsing data from a or string.

    Attributes

    _format     --  Name of the data format that this parses (string, default
                    ""). The format string is a unique identifier for the data
                    format handled by the parser.
    _banks      --  The data from each bank. Each bank contains a (x, y, dx, dy)
                    tuple:
                    x       --  A numpy array containing the independent
                                variable read from the file.
                    y       --  A numpy array containing the profile
                                from the file.
                    dx      --  A numpy array containing the uncertainty in x
                                read from the file. This is None if the
                                uncertainty cannot be read.
                    dy      --  A numpy array containing the uncertainty read
                                from the file. This is None if the uncertainty
                                cannot be read.
    _x          --  Indpendent variable from the chosen bank
    _y          --  Profile from the chosen bank
    _dx         --  Uncertainty in independent variable from the chosen bank
    _dy         --  Uncertainty in profile from the chosen bank
    _meta       --  A dictionary containing metadata read from the file.

    General Metadata

    filename    --  The name of the file from which data was parsed. This key
                    will not exist if data was not read from file.
    nbanks      --  The number of banks parsed.
    bank        --  The chosen bank number.

    """

    _format = ""

    def __init__(self):
        """Initialize the attributes."""
        self._banks = []
        self._meta = {}
        self._x = None
        self._y = None
        self._dx = None
        self._dy = None
        return

    def getFormat(self):
        """Get the format string."""
        return self._format

    def parseString(self, patstring, *args, **kw):
        """Parse a string and set the _x, _y, _dx, _dy and _meta variables.

        When _dx or _dy cannot be obtained in the data format it is set to
        None.

        This wipes out the currently loaded data and selected bank number.
        
        Arguments
        patstring   --  A string containing the pattern

        Raises ParseError if the string cannot be parsed

        """
        raise NotImplementedError()

    def parseFile(self, filename, *args, **kw):
        """Parse a file and set the _x, _y, _dx, _dy and _meta variables.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        filename    --  The name of the file to parse

        Raises IOError if the file cannot be read
        Raises ParseError if the file cannot be parsed

        """
        infile = file(filename, 'r')
        self._banks = []
        self._meta = {}
        filestring = infile.read()
        self.parseString(filestring, *args, **kw)
        infile.close()
        self._meta["filename"] = filename

        if len(self._banks) < 1:
            raise ParseError("There are no data in the banks")

        self.selectBank(0)
        return

    def getNumBanks(self):
        """Get the number of banks read by the parser."""
        return len(self._banks)

    def selectBank(self, index):
        """Select which bank to use.
        
        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end. 

        Arguments:
        index  --  index of bank (integer, starting at 0).

        Raises IndexError if requesting a bank that does not exist

        """
        if index is None: 
            index = self._meta.get("bank", 0)

        numbanks = self.getNumBanks()
        if index > numbanks:
            raise IndexError("Bank index out of range")

        if index < 0:
            index += numbanks

        if index < 0:
            raise IndexError("Bank index out of range")

        self._meta["bank"] = index
        self._meta["nbanks"] = numbanks
        self._x, self._y, self._dx, self._dy = self._banks[index]
        return

    def getData(self, index = None):
        """Get the data.

        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end. 

        Arguments:
        index  --   index of bank (integer, starting at 0). If index is None
                    then the currently selected bank is used.

        This returns (x, y, dx, dy) tuple for the bank. dx is 0 if it cannot
        be determined from the data format.

        """
        self.selectBank(index)

        return self._x, self._y, self._dx, self._dy

    def getMetaData(self):
        """Get the parsed metadata."""
        return self._meta

# End of ProfileParser

class TextParser(ProfileParser):
    """Text parser using numpy.loadtxt.

    Attributes

    _format     --  Name of the data format that this parses ("txt").
    _banks      --  The data from each bank. Each bank contains a (x, y, dx, dy)
                    tuple:
                    x       --  A numpy array containing the independent
                                variable read from the file.
                    y       --  A numpy array containing the profile
                                from the file.
                    dx      --  A numpy array containing the uncertainty in x
                                read from the file. This is None if the
                                uncertainty cannot be read.
                    dy      --  A numpy array containing the uncertainty read
                                from the file. This is None if the uncertainty
                                cannot be read.
    _x          --  Indpendent variable from the chosen bank
    _y          --  Profile from the chosen bank
    _dx         --  Uncertainty in independent variable from the chosen bank
    _dy         --  Uncertainty in profile from the chosen bank
    _meta       --  A dictionary containing metadata read from the file.

    General Metadata

    filename    --  The name of the file from which data was parsed. This key
                    will not exist if data was not read from file.
    nbanks      --  The number of banks parsed.
    bank        --  The chosen bank number.

    Specific Metadata
    The arguments sent to loadtxt, including defaults, are stored in the
    metadata by argument name.

    """

    _format = "txt"

    def parseString(self, patstring, *args, **kw):
        """Use numpy.loadtxt to load data.

        patstring   --  String holding the file contents.

        Arguments are passed to numpy.loadtxt. 
        unpack = True is enforced. 
        The first two arrays returned by numpy.loadtxt are assumed to be x and
        y.  If there is a third array, it is assumed to be dy and if there is a
        fourth it is considered to be dx.  These can be controlled with the
        usecols option. Any other arrays are ignored.

        Raises ParseError if the call to numpy.loadtxt returns fewer than 2
        arrays.

        """
        # Unfortunately, we need to make this into a file-like object again.
        import StringIO
        iofile = StringIO.StringIO(patstring)

        # Enforce unpack
        from diffpy.srfit.util.getcallargs import getcallargs
        import numpy
        callargs = getcallargs(numpy.loadtxt, iofile, *args, **kw)
        callargs["unpack"] = True

        self._meta.update(callargs)

        cols = numpy.loadtxt(**callargs)

        x = y = dy = dx = None
        # Due to using 'unpack', a single column will come out as a single
        # array, thus the second check.
        if len(cols) < 2 or not isinstance(cols[0], numpy.ndarray):
            raise ParseError("numpy.loadtxt returned fewer than 2 arrays")

        x = cols[0]
        y = cols[1]
        if len(cols) > 2:
            dy = cols[2]
        if len(cols) > 3:
            dx = cols[3]

        self._banks.append([x, y, dx, dy])
        return

# End of TextParser

# Register the parser
_registry[TextParser._format] = TextParser


class ParseError(Exception):
    """Exception used by ProfileParsers."""
    pass

__id__ = "$Id$"

# End of ParseError
