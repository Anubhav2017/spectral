/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BSVD.cc.
   Example program that illustrates how to determine the largest 
   singular values of a matrix using arpack++.

   1) Problem description:

      In this example, Arpack++ is called to solve the symmetric problem:

                             (A'*A)*v = sigma*v

      where A is an n by n real band matrix.
      This formulation is appropriate when m >= n.
      The roles of A and A' must be reversed in the case that m < n.

   2) Data structure used to represent the matrix:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 

   3) Included header files:

      File             Contents
      -----------      --------------------------------------------
      bnmatrxw.h       MatrixW, a function that generates matrix A
                       in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arssym.h         The ARSymStdEig class definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arssym.h"
#include "bnmatrxw.h"
#include "arbnsmat.h"
#include "iostream.h"
#include <math.h>


main()
{

  // Defining variables;

  int     i;
  int     n;          // Number of columns in A.
  int     nl;         // Lower bandwidth of A.
  int     nu;         // Upper bandwidth of A.
  double* valA;       // Pointer to an array that stores the elements of A.
  double  cond;       // Condition number of A.
  double* svalue = new double[4];

  // Creating a band matrix with n = 100.

  n  = 100;
  nl = 6;
  nu = 3;
  MatrixW(n, nl, nu, valA);

  // Using ARluNonSymMatrix to store matrix information and to
  // perform the product A'Ax (LU decomposition is not used).

  ARbdNonSymMatrix<double> A(n, nl, nu, valA);

  // Defining what we need: eigenvalues with largest magnitude.

  ARSymStdEig<double, ARbdNonSymMatrix<double> >
    dprob(n, 4L, &A, &ARbdNonSymMatrix<double>::MultMtMv);

  // Finding eigenvalues.

  dprob.Eigenvalues(svalue);

  // Calculating singular values.

  for (i=0; i<4; i++) {
    svalue[i] = sqrt(svalue[i]);
  }

  // Printing some information about the problem.

  cout << endl << "Testing ARPACK++ class ARSymStdEig" << endl;
  cout << "Obtaining singular values by solving (A'*A)*v = sigma*v" << endl;
  cout << endl << "greatest singular values: " << endl; 
  for (i=0; i<4; i++) {
    cout << "  sigma [" << i+1 << "]: " << svalue[i] << endl;
  }

} // main.

