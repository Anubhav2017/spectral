/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE UNSymGSh.cc.
   Example program that illustrates how to solve a real nonsymmetric
   generalized eigenvalue problem in real shift and invert mode using
   the ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      inverse mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.
      The shift sigma is a real number.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The UMFPACK package is called by ARluNonSymGenEig to solve
      some linear systems involving (A-sigma*B). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lnmatrxc.h       StiffnessMatrix, a function that generates
                       matrix A in CSC format.
      lnmatrxd.h       MassMatrix, a function tha generates matrix
                       B in CSC format.
      arunsmat.h       The ARumNonSymMatrix class definition.
      arugnsym.h       The ARluNonSymGenEig class definition.
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

#include "lnmatrxc.h"
#include "lnmatrxd.h"
#include "arunsmat.h"
#include "arugnsym.h"
#include "lnsymsol.h"


main()
{

  // Defining variables;

  int     n;               // Dimension of the problem.
  int     nnza,   nnzb;    // Number of nonzero elements in A and B.
  int     *irowa, *irowb;  // pointers to arrays that store the row
                           // indices of the nonzeros in A and B.
  int     *pcola, *pcolb;  // pointers to arrays of pointers to the
                           // beginning of each column of A and B in
                           // valA and ValB.
  double  rho;             // A parameter used in StiffnessMatrix.
  double  *valA,  *valB;   // pointers to arrays that store the
                           // nonzero elements of A and B.

  // Creating matrices A and B.

  n   = 100;
  rho = 10.0;
  StiffnessMatrix(n, rho, nnza, valA, irowa, pcola);
  ARumNonSymMatrix<double> A(n, nnza, valA, irowa, pcola);

  MassMatrix(n, nnzb, valB, irowb, pcolb);
  ARumNonSymMatrix<double> B(n, nnzb, valB, irowb, pcolb);

  // Defining what we need: the four eigenvectors nearest to 0.0.

  ARluNonSymGenEig<double> dprob(4L, A, B, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

