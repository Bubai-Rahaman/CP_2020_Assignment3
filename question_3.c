/*Fourier transform of sinction using gsl library*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define N 1024
#define PI 3.141592
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

/*Define sinc function*/
float sinc_fun(float x)
{
	if(x==0)
	{
		return(1);
	}
	else
	{
		return(sin(x)/x);
	}
}

int main(void)
{
	float xmin = -500, xmax = 500,dx,pre_fac; //Minimum and maximum value of x is taken as -500 and 500
	dx = (xmax-xmin)/(N-1);
	pre_fac = dx/sqrt(2*PI);
	int i;
	double sampled_data[2*N];
	
	gsl_fft_complex_wavetable * wavetable;
    gsl_fft_complex_workspace * workspace;
    wavetable = gsl_fft_complex_wavetable_alloc (N);
  	workspace = gsl_fft_complex_workspace_alloc (N);
  	
  	float real_factor[N],imag_factor[N],karr[N];
  	
  	/*Sampling data*/
    for (i = 0; i < N; i++)
    {
    	REAL(sampled_data,i) = sinc_fun(xmin+i*dx);
    	IMAG(sampled_data,i) = 0.0;	//The given function is real
    }
	 
	for (i = 0; i < (int) wavetable->nf; i++)
    {
       printf ("# factor %d: %zu\n", i,wavetable->factor[i]);
    }
   
    gsl_fft_complex_forward (sampled_data, 1, N, wavetable, workspace); //Evaluating fft
    
    /*k-values as defined in the documentaion*/
    for (i = 0; i < N; i += 1)
 	{  
   		if(i<(N/2))
   		{
    		karr[i]=(i*2*PI)/(N*dx);
    	}
    	else 
    	{
    		karr[i]=((i-N)*2*PI)/(N*dx);
   		}
    	real_factor[i]=cos(karr[i]*xmin); //xmin is non zero
    	imag_factor[i]=sin(-karr[i]*xmin);
 	}
 	
 	char fname[50];
 	sprintf(fname,"q_3_fftdata.txt"); //output will be stored in q_3_fftdata.txt file
 	FILE *fp=fopen(fname,"w");
 	
 	for (i = 0; i < N; i++)
 	{
 		REAL(sampled_data,i) = pre_fac*(real_factor[i]*REAL(sampled_data,i)-imag_factor[i]*IMAG(sampled_data,i));		//real part
 		IMAG(sampled_data,i) = pre_fac*(real_factor[i]*IMAG(sampled_data,i)+imag_factor[i]*REAL(sampled_data,i));		//imaginary part
 		fprintf(fp,"%f  %f \n",karr[i],REAL(sampled_data,i));
 	}
 	
 	fclose(fp);
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return(0);	
}
