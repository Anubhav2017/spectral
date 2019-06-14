/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LCompGSh.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem in shift and invert mode using
   the ARluCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from a finite element
      discretization of a 1-dimensional convection-diffusion operator
                         (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1], with zero boundary conditions, using
      piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluCompGenEig to solve
      some linear systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxe.h       CompMatrixE, a function that generates matrix
                       A in CSC format.
      lcmatrxf.h       CompMatrixF, a function that generates matrix
                       B in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlgcomp.h       The ARluCompGenEig class definition.
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
#include "lcmatrxe.h"
#include "lcmatrxf.h"
#include "arlnsmat.h"
#include "arlgcomp.h"
#include "lcompsol.h"


main()
{

  // Defining variables;

  int     n;                      // Dimension of the problem.
  int     nnza,   nnzb;           // Number of nonzero elements in A and B.
  int     *irowa, *irowb;         // pointers to arrays that store the row
                                  // indices of the nonzeros in A and B.
  int     *pcola, *pcolb;         // pointers to arrays of pointers to the
                                  // beginning of each column of A and B in
                                  // valA and ValB.
  arcomplex<double> rho;          // parameter used in CompMatrixE.
  arcomplex<double> *valA, *valB; // pointers to arrays that store the
                                  // nonzero elements of A and B.

  // Creating complex matrices A and B.

  n   =  100;
  rho = arcomplex<double>(10.0, 0.0);
  CompMatrixE(n, rho, nnza, valA, irowa, pcola);
  ARluNonSymMatrix<arcomplex<double> > A(n, nnza, valA, irowa, pcola);

  CompMatrixF(n, nnzb, valB, irowb, pcolb);
  ARluNonSymMatrix<arcomplex<double> > B(n, nnzb, valB, irowb, pcolb);

  // Defining what we need: the four eigenvectors nearest to sigma.

  ARluCompGenEig<double> dprob(4L, A, B, arcomplex<double>(1.0,0.0));

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

