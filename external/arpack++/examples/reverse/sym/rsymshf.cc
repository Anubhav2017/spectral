/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE RSymShf.cc.
   Example program that illustrates how to solve a real symmetric
   standard eigenvalue problem in shift and invert mode using the
   ARrcSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and
      invert mode, where A is derived from the central difference
      discretization of the 1-dimensional Laplacian on [0,1] with
      zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      class ARrcSymStdEig requires the user to provide a way to
      perform the matrix-vector product w = OPv, where OP =
      inv[A - sigma*I]. In this example a class called SymMatrixB was
      created with this purpose. SymMatrixB contains a member function,
      MultOPv, that takes a vector v and returns the product OPv in w.

   3) The reverse communication interface:

      This example uses the reverse communication interface, which
      means that the desired eigenvalues cannot be obtained directly
      from an ARPACK++ class.
      Here, the overall process of finding eigenvalues by using the
      Arnoldi method is splitted into two parts. In the first, a
      sequence of calls to a function called TakeStep is combined
      with matrix-vector products in order to find an Arnoldi basis.
      In the second part, an ARPACK++ function like FindEigenvectors
      (or EigenValVectors) is used to extract eigenvalues and
      eigenvectors.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      smatrixb.h       The SymMatrixB class definition.
      arrssym.h        The ARrcSymStdEig class definition.
      rsymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arrssym.h"
#include "smatrixb.h"
#include "rsymsol.h"


template<class T>
void Test(T type)
{

  // Creating a symmetric matrix.

  SymMatrixB<T> B(100,0.0); // n = 100, shift = 0.0.

  // Creating a symmetric eigenvalue problem and defining what we need:
  // the four eigenvectors of B nearest to 0.0.

  ARrcSymStdEig<T> prob(B.ncols(), 4, (T)0.0);

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {

      // Performing matrix-vector multiplication.
      // In shift and invert mode, w = OPv must be performed
      // whenever GetIdo is equal to 1 or -1. GetVector supplies
      // a pointer to the input vector, v, and PutVector a pointer
      // to the output vector, w.

      B.MultOPv(prob.GetVector(), prob.PutVector());

    }
  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvectors();

  // Printing solution.

  Solution(prob);

} // Test.


main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main

