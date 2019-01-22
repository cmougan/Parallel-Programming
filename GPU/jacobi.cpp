#include <cmath>
#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using std::cout;

#ifndef TYPE
#define TYPE double
#endif

#define TOLERANCE 0.001

void
init_simple_diag_dom(int nsize, TYPE* A)
{
  int i, j;

  // In a diagonally-dominant matrix, the diagonal element
  // is greater than the sum of the other elements in the row.
  // Scale the matrix so the sum of the row elements is close to one.

  for (i = 0; i < nsize; ++i) {
    TYPE sum;
    sum = (TYPE)0;
    for (j = 0; j < nsize; ++j) {
      TYPE x;
      x = (rand() % 23) / (TYPE)1000;
      A[i*nsize + j] = x;
      sum += x;
    }
    // Fill diagonal element with the sum
    A[i*nsize + i] += sum;

    // scale the row so the final matrix is almost an identity matrix
    for (j = 0; j < nsize; j++)
      A[i*nsize + j] /= sum;
  }
} // init_simple_diag_dom

int
main(int argc, char **argv)
{
  int nsize; // A[nsize][nsize]
  int i, j, iters, max_iters, riter;
  double start_time, elapsed_time;
  TYPE residual, err;
  TYPE *A, *b, *xnew, *xold;

  // set matrix dimensions and allocate memory for matrices
  nsize = 0;
  if (argc > 1)
    nsize = atoi(argv[1]);
  if (nsize <= 0)
    nsize = 1000;

  max_iters = 0;
  if (argc > 2)
    max_iters = atoi(argv[2]);
  if (max_iters <= 0)
    max_iters = 5000;

  riter = 0;
  if (argc > 3)
    riter = atoi(argv[3]);
  if (riter <= 0)
    riter = 200;

  cout << "nsize = " << nsize << ", max_iters = " << max_iters << "\n";

  A    = new TYPE[nsize*nsize];
  b    = new TYPE[nsize];
  xnew = new TYPE[nsize];
  xold = new TYPE[nsize];

  // generate a diagonally dominant matrix
  init_simple_diag_dom(nsize, A);

  // zero the x vectors, random values to the b vector
  for (i = 0; i < nsize; i++) {
    xnew[i] = (TYPE)0.0;
    xold[i] = (TYPE)0.0;
    b[i] = (TYPE)(rand() % 51) / 100.0;
  }

  start_time = omp_get_wtime();

  // Move diagonal to vector and put zeros
  TYPE *diag = new TYPE[nsize];
  for (i = 0; i < nsize; i++) {
    diag[i] = A[i*nsize+i];
    A[i*nsize+i] = (TYPE)0.0;
  }
  
  //
  // jacobi iterative solver
  //

  residual = TOLERANCE + 1.0;
  iters = 0;
  #pragma acc data copyin(A[:nsize*nsize],b[:nsize],diag[:nsize]) copy(xnew[:nsize],xold[:nsize])
  {
  while ((residual > TOLERANCE*TOLERANCE) && (iters < max_iters)) {
    ++iters;
    if (iters & 1 == 1)
    {
      residual = 0.0;
      #pragma acc kernels loop independent
      for (i = 0; i < nsize; ++i) {
        TYPE rsum = (TYPE)0;
        for (j = 0; j < nsize; ++j) {
          rsum += A[i*nsize + j] * xold[j];
        }
        xnew[i] = (b[i] - rsum) / diag[i];
        TYPE dif = xnew[i] - xold[i];
        residual += dif * dif;
      }
    }
    else {
      residual = 0.0;
      #pragma acc kernels loop independent
      for (i = 0; i < nsize; ++i) {
        TYPE rsum = (TYPE)0;
        for (j = 0; j < nsize; ++j) {
          rsum += A[i*nsize + j] * xnew[j];
        }
        xold[i] = (b[i] - rsum) / diag[i];
        TYPE dif = xnew[i] - xold[i];
        residual += dif * dif;
      }
    }
    if (iters % riter == 0 ) 
      cout << "Iteration " << iters << ", residual is " << sqrt((double)residual) << "\n";
  }}

  elapsed_time = omp_get_wtime() - start_time;
  cout << "\nConverged after " << iters << " iterations and " << elapsed_time << " seconds, residual is " << sqrt((double)residual) << "\n";

  //
  // test answer by multiplying my computed value of x by
  // the input A matrix and comparing the result with the
  // input b vector.
  //
  err = (TYPE)0.0;

  // move diagonal elements
  for (i = 0; i < nsize; i++) {
    A[i*nsize+i] = diag[i];
  }
  delete diag;

  if (iters & 1 == 1)
  {
    for (i = 0; i < nsize; i++) {
      TYPE tmp;
      xold[i] = (TYPE)0.0;
      for (j = 0; j < nsize; j++)
        xold[i] += A[i*nsize + j] * xnew[j];
      tmp = xold[i] - b[i];
      err += tmp * tmp;
    }
  } else
  {
    for (i = 0; i < nsize; i++) {
      TYPE tmp;
      xnew[i] = (TYPE)0.0;
      for (j = 0; j < nsize; j++)
        xnew[i] += A[i*nsize + j] * xold[j];
      tmp = xnew[i] - b[i];
      err += tmp * tmp;
    }
  }
  err = sqrt((double)err);
  cout << "Solution error is " << err << "\n";
  if (err > TOLERANCE)
    cout << "****** Final Solution Out of Tolerance ******\n" << err << " > " << TOLERANCE << "\n";

  delete A;
  delete b;
  delete xnew;
  delete xold;
  return 0;
}

