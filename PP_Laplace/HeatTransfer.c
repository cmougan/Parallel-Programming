#include<stdio.h>
#include<math.h>


double T,time;
int N;


int main(int argc, char *argv[])
{
	printf("Indroduce the T (Temperature) and N (dimension): ");
    scanf("%lf",&T);
 	scanf("%d",&N);
 	printf("Indroduce the t (time when the measure is wanted): ");
    scanf("%lf",&time);
 	
	int i,x,t;
    double alpha=1.0,L=0.345678;
	double **U=malloc(N * sizeof(double));
 	for (i=0;i<N;++i)
	{
    	U[i] = malloc(2 * sizeof(double));
	}
	double difx=L/N,dift=0.4*difx*difx; //(dift/difx^2)<=0.5
	int time_cicles=(int)ceil(time/dift);
	
	for(x=0;x<N;x++){
        U[x][0]=0.0;
		U[x][1]=0.0;
    }
    U[0][0]=T;
    U[N-1][0]=T;
    U[0][1]=T;
    U[N-1][1]=T;
    
    for(x=0;x<N;x++){
        printf(" (%e,%e) ",U[x][0],U[x][1]);
    }
    int iter;
    for(iter=0;iter<time_cicles;iter++){		
		x=0;
		do{
			x++;
			U[x][1]=alpha*dift*((U[x+1][0]-2.0*U[x][0]+U[x-1][0])/(difx*difx))+U[x][0];
		}while(x<(N-2));
		
		for(x=0;x<N;x++){
	            	U[x][0]=U[x][1];	
		}
	}
	
	printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
	for(x=0;x<N;x++){
                printf(" (%e,%e) ",U[x][0],U[x][1]);
    }
	return 0;
}
