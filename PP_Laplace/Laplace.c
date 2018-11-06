#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int  n, m;
double tol=0.001 ;
double **A, **Anew; 																			// Al hacerlas globales estan inicializadas a cero por defecto.

void main(int argc, char**argv)
{
   /*
    * Metodo con un argumento de entrada que es int_max en caso de que no haya argumento de entrada default=200 
    * y se piden n y m por pantalla
    * 
    */
  int iter_max=(argv[1]!=0)?( atoi(argv[1])):(int)200;
  printf("Indroduce the n and m dimensions: ");
  scanf("%d",&n);
  scanf("%d",&m);
  
  
  int i, j,converged=0,iter_num=0;																// Se declar i j y se inicializa converge e iter_num a cero
  double error,error_aux ;
  A=malloc(n * sizeof(double));
  Anew=malloc(n * sizeof(double));																// Malloc deja decidir al usuario, que notifica cuantas columnas hay
  for (i=0;i<n;++i)																				// Malloc deja decidir al usuario, se notifica cuantas filas hay
	{
    	A[i] = malloc(m * sizeof(double));
    	Anew[i] = malloc(m * sizeof(double));
	}
	/*
	 * Condiciones de Contorno
	 * mediante un unico bucle
	 */	
  for(i=0;i<n;i++){
	    A[i][0]= sin(i*M_PI/(m-1));
	    Anew[i][0]= A[i][0];
	    A[i][m-1]=sin(i*M_PI/(m-1))*exp(-M_PI);
	    Anew[i][m-1]= A[i][m-1];
  }
  /*
   * Calculo de diferencias finitas entre dos matrices Anew es la iteracion calculada
   * y tras cada ciclo se igual la post a la actual.
   */
  while(!converged && iter_num<iter_max){
    error=0.0;
    iter_num++;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            if(!(i==0||i==n-1||j==0||j==m-1)){
	            Anew[i][j]=(A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1])/4.0;
	            error_aux=Anew[i][j]-A[i][j];
	            error=(error_aux>error)?error_aux:error;
	        }
        }
    }
    for(i=0;i<n;i++){
	  for (j=0;j<m;j++){
	    A[i][j]=Anew[i][j];
 	  }
    } 
    //Printeo por pantalla de la iteracion y el error existente de la misma
    if(iter_num%10==0){
        printf("In the iter %i , te error is: %e \n",iter_num,error);
    }
    //Comprobación de del error para salir del ciclo while.
    converged=((error=sqrt(error))<tol)?1:0;
  }
  
  //Printeo de los resultados.
  printf(" The result is: \n ");
  for(i=0;i<n;i++){
	for (j=0;j<m;j++){
	  printf(" %lf ",A[i][j]);
	}
	printf("\n ");
  }
  printf("obtained in the iter: %i . With max error: %e \n ",iter_num, error);
}
