/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE HNSymGen.cc.
   Example program that illustrates how to solve a real
   nonsymmetric generalized eigenvalue problem using the
   ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in
      regular or shift and invert mode, where A is read
      from a file in Harwell-Boeing format.

   2) Included header files:

      File             Contents
      -----------      -------------------------------------------
      arerror.h        The ArpackError class definition.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlgnsym.h       The ARluNonSymGenEig class definition.
      lnsymsol.h       The Solution function.

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
#include "arlnsmat.h"
#include "arlgnsym.h"
#include "lnsymsol.h"


void PrintHelp()
/*
  Prints information about hnsymgen usage.
*/
{

  cout << "ARPACK++ version 1.0 may 97" << endl;
  cout << "usage:   hnsymgen [parameters] file" << endl;
  cout << "parameters:" << endl;
  cout << "      -n <number of desired eigenvalues>" << endl;
  cout << "      -c <number of Arnoldi vectors per iteration>" << endl;
  cout << "      -l <maximum number of iterations>" << endl;
  cout << "      -s <real part of the shift>" << endl;
  cout << "      -i <imaginary part of the shift>" << endl;
  cout << "      -p <part of inv(A-sB)*v considered (R or I))>" << endl;
  cout << "      -t <stopping criterion>" << endl;
  cout << "      -u <LU pivot threshold>" << endl;
  cout << "      -o <column ordering for factorization>" << endl;
  cout << "      -w <desired portion of the spectrum. " << endl;
  cout << "          acceptable values: LM, SM, LR, SR, LI, SI>" << endl;
  cout << endl;

} // PrintHelp.


bool ReadParameters(int n, char* v[], int &nev, int &ncv, int &maxit,
                    int &order, bool &shift, double &sigmar, 
                    double &sigmai, char &part, double &tol, double &thresh, 
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
  order   = 1;
  shift   = false;
  sigmar  = 0.0;
  sigmai  = 0.0;
  tol     = 0.0;
  thresh  = 0.1;
  part    = 'R';
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
          sigmar = atof(v[i++]);
          shift  = true;
          break;
        case 'i':
          sigmai = atof(v[i++]);
          shift  = true;
          break;
        case 'p':
          part = v[i++][0];
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
  double sigmar;
  double sigmai;
  double tol;
  double thresh;
  char   part;
  char*  which;
  char*  fileA;
  char*  fileB;

  // Reading parameters.

  if (!ReadParameters(argc, argv, nev, ncv, maxit, order, shift, sigmar, 
                      sigmai, part, tol, thresh, which, fileA, fileB)) {
    return 0;
  }

  // Reading and storing matrices A and B.

  ARluNonSymMatrix<double> A(fileA, thresh, order);
  ARluNonSymMatrix<double> B(fileB, thresh, order);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<double> dprob(nev, A, B, which, ncv, tol, maxit);

  // Defining the shift.

  if (shift) {
    if (sigmai == 0.0) {
      dprob.SetShiftInvertMode(sigmar);
    }
    else {
      dprob.SetComplexShiftMode(part, sigmar, sigmai);
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

