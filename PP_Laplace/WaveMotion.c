#include<stdio.h>
#include<math.h>

int  X,T;


int main(int argc, char *argv[])
{
	printf("Indroduce  X (dimension): ");
    scanf("%d",&X);
 	printf("Indroduce the T (time intervals): ");
    scanf("%d",&T);
 	
	int i,x;
    double L=0.345678;
    
	double **U=malloc((X+1)* sizeof(double));
 	for (i=0;i<(X+1);i++)
	{
    	U[i] = malloc((T+2) * sizeof(double));
	}
	for(x=0;x<X+1;x++){
        U[x][0]=sin(x*M_PI/X);
		U[x][1]=sin(x*M_PI/X)*cos(M_PI/T);
    }
    
    for(i=1;i<(T+1);i++){
		U[0][i+1]=0;
		U[X][i+1]=0;	
		x=0;
		do{
			x++;
			U[x][i+1]=L*(U[x-1][i]-2*U[x][i]+U[x+1][i])+2*U[x][i]-U[x][i-1];
		}while(x<(X-1));
		
	}
	for(i=0;i<(T+2);i++){
		for(x=0;x<(X+1);x++){
			printf(" %f \t",U[x][i]);
		}
		printf("\n");
	}
}
