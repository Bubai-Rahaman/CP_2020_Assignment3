/*Fourier transform of Gaussian function using FFTW library*/

#include<stdlib.h>
#include<stdio.h>
#include <fftw3.h>
#include<math.h>
#define N 1024
#define PI 3.14159265

/*Define Gaussian function*/
float Gaussian(float x)
{ 
	return(exp(-x*x));
}

void main()
{
 
  int i;
  float xmin=-50.0,xmax=50.0,dx; //setting the x range
  dx=(xmax-xmin)/(N-1);
  float pre_fac = (dx/sqrt(2*PI));
  
  fftw_complex *sampled_data, *fft_data ,*factor;
  fftw_plan p;
  
  sampled_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fft_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  factor=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  
  float xarr[N], karr[N];
  
  /*sampling the data*/
  for (i = 0; i < N; i++)
  {  
    xarr[i]=xmin+i*dx;
    sampled_data[i][0]=Gaussian(xarr[i]);
    sampled_data[i][1]=0.0; //Gaussian is real
  }
  
 /*calculating dft*/
  p = fftw_plan_dft_1d(N, sampled_data, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  
  /*k-values*/
  for (i = 0; i < N; i++)
  {  
   if(i<(N/2))
   {
 		karr[i]=(i/(N*dx))*(2*PI);
 	}
    else
    {
    	karr[i]=((i-N)/(N*dx))*(2*PI);
    }
   
    factor[i][0]=cos(karr[i]*xmin); 	//xmin is non zero
    factor[i][1]=sin(-karr[i]*xmin);
  }
 
 char fname[50];
 sprintf(fname,"q_4_fftdata.txt"); // output will be stored in q_4_fftdata.txt
 FILE *fp=fopen(fname,"w");
  for (i = 0; i < N; i++)
  {  
   fft_data[i][0] = pre_fac*(factor[i][0]*fft_data[i][0]-factor[i][1]*fft_data[i][1]); //real part
   fft_data[i][1] = pre_fac*(factor[i][0]*fft_data[i][1]+factor[i][1]*fft_data[i][0]);	//imaginary part
   
   fprintf(fp,"%f  %f \n",karr[i],fft_data[i][0]);
  }
  
  fclose(fp);
  fftw_destroy_plan(p);
  fftw_free(sampled_data);
  fftw_free(fft_data);
  fftw_free(factor);
}
