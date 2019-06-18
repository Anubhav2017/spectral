# ARPACK++ v1.0 8/1/1997
# c++ interface to ARPACK code.
# This file contains some definitions used to compile arpack++ examples
# with the g++ compiler on a SUN workstation.


# Defining the machine.

MACHINE      = SUN4

# Defining the compiler.

CPP          = g++

# Defining ARPACK++ directories.
# ARPACKPP_INC is the directory that contains all arpack++ header files.
# SUPERLU_DIR and UMFPACK_DIR must be set to ARPACKPP_INC.

ARPACKPP_DIR = $(HOME)/arpack++
ARPACKPP_INC = $(ARPACKPP_DIR)/include
SUPERLU_DIR  = $(ARPACKPP_INC)
UMFPACK_DIR  = $(ARPACKPP_INC)

# Defining ARPACK, LAPACK, UMFPACK, SUPERLU, BLAS and FORTRAN libraries.
# See the arpack++ manual or the README file for directions on how to 
# obtain arpack, umfpack and SuperLU packages. 
# UMFPACK_LIB and SUPERLU_LIB must be declared only if umfpack and superlu 
# are going to be used. BLAS and LAPACK fortran routines are distributed 
# along with arpack fortran code, but the user should verify if optimized 
# versions of these libraries are available before installing arpack. 
# The fortran libraries described below are those required to link
# fortran and c++ code together in a Sun SparcStation. Other libraries
# should be defined if the user intends to compile arpack++ on another
# environment.

ARPACK_LIB   = $(HOME)/ARPACK/libarpack.a
LAPACK_LIB   = $(HOME)/LAPACK/liblapack.a
UMFPACK_LIB  = $(HOME)/UMFPACK/libumfpack.a
SUPERLU_LIB  = -lm $(HOME)/SuperLU/superlu.a
BLAS_LIB     = $(HOME)/LAPACK/libblas.a
FORTRAN_LIBS = -lm /opt/SUNWspro/SC4.2/lib/libF77.a \
               -lm /opt/SUNWspro/SC4.2/lib/libsunmath.a \
               -lm /opt/SUNWspro/SC4.2/lib/libM77.a \
               -lm /opt/SUNWspro/SC4.2/lib/libsunmath.a 

# Defining g++ flags.

CPP_WARNINGS = -Wall
CPP_DEBUG    = -g -fhandle-exceptions
CPP_OPTIM    = -O -fhandle-exceptions
CPP_LIBS     =

CPP_FLAGS    = $(CPP_DEBUG) -D$(MACHINE) -I$(ARPACKPP_INC) 

# Putting all libraries together.

ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) $(LAPACK_LIB) \
               $(BLAS_LIB) $(FORTRAN_LIBS)

# defining paths.

vpath %.h  $(ARPACK_INC)
