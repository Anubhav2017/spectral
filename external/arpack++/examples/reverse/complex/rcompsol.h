/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE RCompSol.h
   Printing eigenvalues of a complex problem
   (ARrcCompStdEig version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RCOMPSOL_H
#define RCOMPSOL_H

#include "arrscomp.h"

template<class FLOAT>
void Solution(ARrcCompStdEig<FLOAT> &Prob)
/*
  Prints eigenvalues on standard "cout" stream.
*/

{

  int i, n, nconv, mode;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARrcCompStdEig " << endl;
  cout << "Complex eigenvalue problem: A*x - lambda*x" << endl;
  switch (mode) {
  case 1:
    cout << "Regular mode" << endl << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl << endl;
  }

  cout << "Dimension of the system            : " << n              << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
  cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {
    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << endl;
    }
    cout << endl;
  }

} // Solution


#endif // RCOMPSOL_H
