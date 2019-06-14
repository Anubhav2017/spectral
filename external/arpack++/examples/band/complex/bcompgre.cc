/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BCompGRe.cc.
   Example program that illustrates how to solve a complex generalized
   eigenvalue problem in regular mode using the ARluCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular mode,
      where A and B are derived from the finite element discretization
      of the 1-dimensional convection-diffusion operator
                       (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1], with zero boundary conditions, using
      piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 
      {ndiagL, ndiagU, B}: matrix B in band format.

   3) Library called by this example:

      The LAPACK package is called by ARluCompGenEig to solve
      some linear systems involving B.

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      bcmatrxb.h       CompMatrixE, a function that generates matrix
                       A in band format.
      bcmatrxc.h       CompMatrixF, a function that generates matrix
                       B in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arbgcomp.h       The ARluCompGenEig class definition.
      lcompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "bcmatrxb.h"
#include "bcmatrxc.h"
#include "arbnsmat.h"
#include "arbgcomp.h"
#include "lcompsol.h"


main()
{

  // Defining variables;

  int               n;            // Dimension of the problem.
  int               ndiagL;       // Lower bandwidth of A and B.
  int               ndiagU;       // Upper bandwidth of A and B.
  arcomplex<double> rho;          // Parameter used to define A.
  arcomplex<double> *valA, *valB; // pointers to arrays that store
                                  // the elements of A and B.

  // Creating complex matrices A and B.

  n   =  100;
  rho = arcomplex<double>(10.0, 0.0);
  CompMatrixB(n, rho, ndiagL, ndiagU, valA);
  ARbdNonSymMatrix<arcomplex<double> > A(n, ndiagL, ndiagU, valA);

  CompMatrixC(n, ndiagL, ndiagU, valB);
  ARbdNonSymMatrix<arcomplex<double> > B(n, ndiagL, ndiagU, valB);

 // Defining what we need: the four eigenvectors with largest magnitude.

  ARluCompGenEig<double> dprob(4L, A, B);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.
