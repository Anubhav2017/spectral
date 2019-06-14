/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BNSymGRe.cc.
   Example program that illustrates how to solve a real nonsymmetric
   generalized eigenvalue problem in regular mode using the
   ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular
      mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.

   2) Data structure used to represent matrices A and B:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 
      {ndiagL, ndiagU, B}: matrix B in band format.

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymGenEig to solve
      some linear systems involving B.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      bnmatrxb.h       StiffnessMatrix, a function that generates
                       matrix A in band format.
      bnmatrxc.h       MassMatrix, a function that generates 
                       matrix B in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arbgnsym.h       The ARluNonSymGenEig class definition.
      lnsymsol.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "bnmatrxb.h"
#include "bnmatrxc.h"
#include "arbnsmat.h"
#include "arbgnsym.h"
#include "lnsymsol.h"


main()
{

  // Defining variables;

  int     n;       // Dimension of the problem.
  int     ndiagL;  // Lower bandwidth of A and B.
  int     ndiagU;  // Upper bandwidth of A and B.
  double  rho;     // Parameter used to define A.
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  // Creating matrices A and B.

  n   = 100;
  rho = 10.0;
  StiffnessMatrix(n, rho, ndiagL, ndiagU, valA);
  ARbdNonSymMatrix<double> A(n, ndiagL, ndiagU, valA);

  MassMatrix(n, ndiagL, ndiagU, valB);
  ARbdNonSymMatrix<double> B(n, ndiagL, ndiagU, valB);

  // Defining what we need: the four eigenvectors with largest magnitude.

  ARluNonSymGenEig<double> dprob(4L, A, B);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

