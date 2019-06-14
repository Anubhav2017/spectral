/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE HSymGen.cc.
   Example program that illustrates how to solve a real
   symmetric generalized eigenvalue problem using the
   ARluNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in
      regular, shift and invert buckling or Cayley mode,
      where A is read from a file in Harwell-Boeing format.

   2) Included header files:

      File             Contents
      -----------      -------------------------------------------
      arerror.h        The ArpackError class definition.
      arlsmat.h        The ARluSymMatrix class definition.
      arlgsym.h        The ARluSymGenEig class definition.
      lsymsol.h        The Solution function.

   3) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include <iostream.h>
#include "arerror.h"
#include "arlsmat.h"
#include "arlgsym.h"
#include "lsymsol.h"


void PrintHelp()
/*
  Prints information about hsymgen usage.
*/
{

  cout << "ARPACK++ version 1.0 may 97" << endl;
  cout << "usage:   hsymgen [parameters] file" << endl;
  cout << "parameters:" << endl;
  cout << "      -n <number of desired eigenvalues>" << endl;
  cout << "      -c <number of Arnoldi vectors per iteration>" << endl;
  cout << "      -l <maximum number of iterations>" << endl;
  cout << "      -s <shift>" << endl;
  cout << "      -i <invert mode: 'S'hift-invert, 'B'uckling or 'C'ayley)>";
  cout << endl;
  cout << "      -t <stopping criterion>" << endl;
  cout << "      -u <LU pivot threshold>" << endl;
  cout << "      -o <column ordering for factorization>" << endl;
  cout << "      -w <desired portion of the spectrum. " << endl;
  cout << "          acceptable values: LM, SM, LA, SA, BE>" << endl;
  cout << endl;

} // PrintHelp.


bool ReadParameters(int n, char* v[], int &nev, int &ncv, int &maxit,
                    int &order, bool &shift, double &sigma,
                    char &mode, double &tol, double &thresh,
                    char* &which, char* &fileA, char* &fileB)
/*
  Reads parameters from the command line.
*/
{

  int  i;
  bool ok;

  // Defining default values.

  nev     = 5;
  ncv     = 0;
  maxit   = 0;
  order   = 2;
  shift   = false;
  sigma   = 0.0;
  tol     = 0.0;
  thresh  = 0.1;
  mode    = 'S';
  which   = "LM";
  fileA   = " ";
  fileB   = " ";
  ok      = true;

  // Returning if the number of parameters is even.

  if ((n==1)||(!(n%2))) {
    ok = false;
  }
  else {

    // Reading parameters.

    i = 1;
    while ((i<(n-2)) && (ok)) {
    
      if ((v[i][0] != '-') || (strlen(v[i]) != 2)) {
        ok = false;
        i += 2;
      }
      else {
        switch (v[i++][1]) {
        case 'n':
          nev = atoi(v[i++]);
          break;
        case 's':
          sigma = atof(v[i++]);
          shift = true;
          break;
        case 'i':
          mode = v[i++][0];
          if ((mode!='S')&&(mode!='B')&&(mode!='C')) {
            cout << "invalid invert mode: " << mode << endl;
            ok = false;
          }
          break;
        case 'w':
          which = v[i++];
          break;
        case 'c':
          ncv = atoi(v[i++]);
          break;
        case 't':
          tol = atof(v[i++]);
          break;
        case 'l':
          maxit = atoi(v[i++]);
          break;
        case 'o':
          order = atoi(v[i++]);
          break;
        case 'u':
          thresh = atof(v[i++]);
          break;
        default :
          cout << "unrecognized parameter: -" << v[i-1][1] << endl;
          ok = false;
        }
      }
    }
  }

  if (ok) {
    fileA = v[i++];
    fileB = v[i];
  }
  else {
    PrintHelp();
  }

  return ok;

} // ReadParameters.


main(int argc, char* argv[])
{

  // Defining variables.

  int    nev;
  int    ncv;
  int    maxit;
  int    order;
  bool   shift;
  double sigma;
  double tol;
  double thresh;
  char   mode;
  char*  which;
  char*  fileA;
  char*  fileB;

  // Reading parameters.

  if (!ReadParameters(argc, argv, nev, ncv, maxit, order, shift, sigma,
                      mode, tol, thresh, which, fileA, fileB)) {
    return 0;
  }

  // Reading and storing matrices A and B.

  ARluSymMatrix<double> A(fileA, thresh, order);
  ARluSymMatrix<double> B(fileB, thresh, order);

  // Defining the eigenvalue problem.

  ARluSymGenEig<double> dprob(nev, A, B, which, ncv, tol, maxit);

  // Defining the shift.

  if (shift) {
    switch (mode) {
    case 'S':
      dprob.SetShiftInvertMode(sigma);
      break;
    case 'B':
      dprob.SetBucklingMode(sigma);
      break;
    case 'C':
      dprob.SetCayleyMode(sigma);
    }
  }

  // Finding eigenvalues and eigenvectors.

  try {
    dprob.FindEigenvectors();
  }
  catch (ArpackError) { return 0; }

  // Printing solution.

  Solution(A, B, dprob);

} // main.

