# figure out user variables
#
ifdef DEBUG
CFLAGS = $(DEBUG) -fPIC -Wall
else
CFLAGS = -fPIC -O3 -w -ffast-math -pipe -funroll-loops
endif

# PYTHON_VERSION
ifndef PYTHON_VERSION
PYTHON_VERSION = $(shell python -c 'import sys; print sys.version[:3]')
endif

# PYTHON_INCDIR
ifndef PYTHON_INCDIR
PYTHON_INCDIR = /usr/include/python$(PYTHON_VERSION)
endif

# NUMPY_INCDIR
ifndef NUMPY_INCDIR
NUMPY_INCDIR = $(shell python -c "import numpy; import os.path; print os.path.split(numpy.__file__)[0] + '/core/include'")
endif

CC        = g++
LINK      = g++
LINKFLAGS = -L$(LIB) -L.

OBJ = Visitor.o Callout.o Evaluator.o Literal.o Operator.o UFuncOperator.o Argument.o Clicker.o

lib : $(OBJ)
	$(LINK) -shared $(LINKFLAGS) $^ -o libsrfit.so

main : main.o $(OBJ)
	$(CC) $^ -o $@ -lboost_python -lpython$(PYTHON_VERSION) 

%.o:%.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(PYTHON_INCDIR) -I$(NUMPY_INCDIR)

.PHONY : clean
clean: 
	-rm -f *.pyc *.so *.o *.a
