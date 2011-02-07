"""Utilities for testing."""

class TestViewer(object):
    def __init__(self):
        self.msg = None
    def _respond(self, msg):
        self.msg = msg

# end class TestViewer

class TestVisitor(object):

    def __init__(self):
        self.nodetype = None

    def onParameter(self, obj):
        self.nodetype = "parameter"

    def onOperator(self, obj):
        self.nodetype = "operator"

    def onContainer(self, obj):
        self.nodetype = "container"

# end class TestVisitor

class TestListContainer(list):

    def __init__(self):
        self.a = 1
        self.b = 2
        self._c = 3
        self.extend([4,5])
        return

    def configure(self):
        self.configured = True

    def getc(self):
        return self._c

    def getpar(self, name):
        return getattr(self, name)

    def getwhatever(self, ignored):
        return self.a

    def calc(self, d):
        return self.a + self.b + self._c + d

# End class TestListContainer

class TestDictContainer(dict):

    def __init__(self):
        self.a = 1
        self.b = 2
        self._c = 3
        self["d"] = 4
        self["e"] = 5
        return

    def configure(self):
        self.configured = True

    def getc(self):
        return self._c

    def getpar(self, name):
        return getattr(self, name)

    def calc(self, d):
        return self.a + self.b + self._c + d

# End class TestDictContainer
