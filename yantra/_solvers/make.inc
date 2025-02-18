ifndef COMPILER
COMPILER = GCC
endif

ifeq ($(COMPILER),GCC)
FC = GNU95
F90FLAGS = "-fopenmp" -lgomp 
OPTFLAGS ="-O3"
XF2PYFLAGS1 +="--quiet"
XF2PYFLAGS2 +="--quiet"
endif

ifeq ($(COMPILER),Intel)
FC = intelem
F90FLAGS = "-openmp" -lguide 
OPTFLAGS ="-O3"
XF2PYFLAGS1 +="--quiet"
XF2PYFLAGS2 +="--quiet"
endif
