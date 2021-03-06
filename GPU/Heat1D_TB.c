#include <stdio.h>
#include <stdlib.h>

// Define here as constant for easy change
#define REAL float

void printV (REAL A[], int N)
{
  int x;
  for (x=0; x<=N+1; x++)
  {
    printf("%d: %1.10f\n", x, A[x]);
  }
}

int main(int argc, char **argv)
{
  int  x, t, N= 10000000, T=100;
  REAL L= 0.001234, L2, S,Lft, Md, Rgt;;
  REAL *U1, *U2;

  if (argc>1) { T = atoi(argv[1]); } // get  first command line parameter
  if (argc>2) { N = atoi(argv[2]); } // get second command line parameter
  if (argc>3) { L = atof(argv[3]); } // get  third command line parameter
 
  if (N < 1 || T < 1 || L >= 0.5) {
    printf("arguments: T N L (T: steps, N: vector size, L < 0.5)\n");
    return 1;
  }

  U1 = (REAL *) malloc ( sizeof(REAL)*(N+2) );
  U2 = (REAL *) malloc ( sizeof(REAL)*(N+2) );
  if (!U1 || !U2) { printf("Cannot allocate vectors\n"); exit(1); }

  // initialize temperatures at time t=0  
  for (x=0; x<=N+1; x++)
    U1[x] = U2[x] = 0.0;
 
  // initialize fixed boundary conditions on U1 and U2
  {
    U1[0]  = U2[0]  = 1.0;
    U1[N+1]= U2[N+1]= 9.0;
  }

  // Set initial temperature on two special points
  {
    U1[N/3]  = 8.0;
    U1[4*N/7]= 3.0;
  }

  printf("Stencil computation of %d steps on 1-D vector of %d elements with L=%1.10e\n", T, N, L);

  L2 = (1.0f-2.0f*L);
  #pragma acc data copyin(U1[:(N+2)],U2[:(N+2)])
  {
  for (t=1; t<=T/2; t++)  // loop on time
    if ((t&1) == 1) {// t is odd
	  // prologue: special case for x==1
	  Lft =    U1[0];
	  Md  =    L2*U1[1] + L*(U1[2]+U1[0]);
	  Rgt =    L2*U1[2] + L*(U1[3]+U1[1]);
	  U2[1] = L2*Md + L*(Rgt+Lft);
	  // Normal cases
  #pragma acc kernels loop independent
	  for (x=2; x<N; x++) {
	    Lft =    L2*U1[x-1] + L*( U1[x] +U1[x-2]);
	    Md  =    L2*U1[x]   + L*(U1[x+1]+U1[x-1]);
	    Rgt =    L2*U1[x+1] + L*(U1[x+2]+ U1[x]);
	    U2[x] = L2*Md + L*(Rgt+Lft);
	  }
	  // epilogue: special case for x==N
	  Lft =    L2*U1[N-1] + L*( U1[N] +U1[N-2]);
	  Md  =    L2*U1[N]   + L*(U1[N+1]+U1[N-1]);
	  Rgt =    U1[N+1];
	  U2[N] = L2*Md + L*(Rgt+Lft);
	 }
   else{            // t is even
    // prologue: special case for x==1
	  Lft =    U2[0];
	  Md  =    L2*U2[1] + L*(U2[2]+U2[0]);
	  Rgt =    L2*U2[2] + L*(U2[3]+U2[1]);
	  U1[1] = L2*Md + L*(Rgt+Lft);
	  // Normal cases
   #pragma acc kernels loop independent
	  for (x=2; x<N; x++) {
	    Lft =    L2*U2[x-1] + L*( U2[x] +U2[x-2]);
	    Md  =    L2*U2[x]   + L*(U2[x+1]+U2[x-1]);
	    Rgt =    L2*U2[x+1] + L*(U2[x+2]+ U2[x]);
	    U1[x] = L2*Md + L*(Rgt+Lft);
	  }
	  // epilogue: special case for x==N
	  Lft =    L2*U2[N-1] + L*( U2[N] +U2[N-2]);
	  Md  =    L2*U2[N]   + L*(U2[N+1]+U2[N-1]);
	  Rgt =    U2[N+1];
	  U1[N] = L2*Md + L*(Rgt+Lft);
	 }

  if ((t&1) == 1)
  {
    if (N<=100)
      printV(U1,N);
    else
    {S=0.0;
 #pragma acc kernels loop independent
      for ( x=0; x<=N+1; x++) // compute checksum of final state 
         S = S + U1[x];
      printf("\nCheckSum = %1.10e\n", S); 
     }
  } 
  else
  {
    if (N<=100)
      printV(U2,N);
    else
    {S=0.0;
 #pragma acc kernels loop independent
      for ( x=0; x<=N+1; x++) // compute checksum of final state 
         S = S + U2[x];
      printf("\nCheckSum = %1.10e\n", S);
     }
  }
 }
}

