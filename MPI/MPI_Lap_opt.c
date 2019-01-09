#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define BUFSIZE 20000

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n, int m, int sw_inner)
{
  int i, j;
  float error=0.0f;
  if (sw_inner == 1){

  for ( j=2; j < m-2; j++ )
    for ( i=1; i < n-1; i++ )
    {
      out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
      error = max_error( error, out[j*n+i], in[j*n+i] );
    }
  }
  
  if (sw_inner == 0){
  
  for ( i=1; i < n-1; i++ )
    {
      out[1*n+i]= stencil(in[1*n+i+1], in[1*n+i-1], in[(1-1)*n+i], in[(1+1)*n+i]);
      out[(m-1)*n+i]= stencil(in[(m-1)*n+i+1], in[(m-1)*n+i-1], in[((m-1)-1)*n+i], in[((m-1)+1)*n+i]);
      error = max_error( error, out[(m-1)*n+i], in[(m-1)*n+i] );
      error = max_error( error, out[1*n+i], in[1*n+i] );
    }

}

  return error;
}


void laplace_init(float *in, int n, int m)  //Same problem as before
{
  int i,proc;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in, 0, (n*(m+2))*sizeof(float));
  for (i=0; i<n; i++) {
    float V = in[i*m] = sinf(pi*i / (n-1));
    in[ i*m+n-1 ] = V*expf(-pi);
  }
}

int main(int argc, char** argv)
{
  int n = 512,rank,size, m,proc;
  int iter_max = 100;
  float *A, *temp, *Aprim, *Aprimlast;
    
  const float tol = 1.0e-3f;
  float error= 1.0f, global_error = 1.00f, error_step;
  float new_error = 1.0f;
  float t1,t2, diftime;    

 
  

  MPI_Init(&argc, &argv);
   
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  
   
  MPI_Status status;
  MPI_Request request;
  t1 = MPI_Wtime();
 
  // get runtime arguments 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }
   
  
  m = n/size; 
  printf("%d\n",m);
	
  A    = (float*) malloc( (n*(m+2))*sizeof(float) ); // We set up a matrix of dim n*(m+2) where m is n/total number of processor, and +2 is for the auxilliary space 
  temp = (float*) malloc( (n*(m+2))*sizeof(float) ); 
  Aprim = (float*) malloc(n*sizeof(float));
  Aprimlast = (float*) malloc(n*sizeof(float));
  //  set boundary conditions
  laplace_init (A, n, m);
  laplace_init (temp, n, m);
  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations, and  %d rank\n", 
         n, n, iter_max,rank );

   


   
  
  int iter = 0;
  int i = 0;

  while ( global_error > tol*tol && iter < iter_max ) //Modification to optimize the waiting time in the recieves 
  { 
    iter++;
    
        memcpy(Aprim,&A[n],sizeof(n));
	memcpy(Aprimlast,&A[(n*(m-1))+1],sizeof(n));
   
    if (rank > 0) {
		MPI_Isend(&Aprim[0],n,MPI_FLOAT,rank-1,1, MPI_COMM_WORLD, &request); //We send the "second" row of the matrix 
		} 
    if (rank < size - 1 ) {
		MPI_Isend(&Aprimlast[0],n,MPI_FLOAT,rank+1,2, MPI_COMM_WORLD, &request); //We send the "penultima" row of the matrix 
		}
	error= laplace_step (A, temp, n, m,1);
   if (rank > 0) {
		MPI_Recv(&A[0],n,MPI_FLOAT,rank-1,2, MPI_COMM_WORLD,&status);  //We receive the "first" (auxilliary) row of the matrix	
		} 
    if (rank < size - 1 ) {
		MPI_Recv(&A[(n*m)+1],n,MPI_FLOAT,rank+1,1, MPI_COMM_WORLD,&status); //We receive the last (auxiliary) row of the matrix
		}

	
    error_step = laplace_step (A, temp, n, m,0);
    error>error_step? error_step:error;
    
    MPI_Reduce(&error,&global_error,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD); // Here, the error computed by each node is send to $global_error of node 0, once this node reaches this node, it skips from the while  
    
		
	
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
  }
  
  
  free(A); free(temp);
  t2 = MPI_Wtime();
  MPI_Finalize();
  
  if (rank == 0){
  diftime = t2-t1;
  error = sqrtf(global_error);
  printf("Total Iterations: %5d, ERROR: %0.6f and rank %d (%0.6f seconds elapsed)\n", iter, error, rank, diftime);}

}
