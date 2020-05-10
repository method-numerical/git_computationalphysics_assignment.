#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h> 	 	

#include <gsl/gsl_fft_complex.h>

void main()
{
int n=512;
double dx=1;
double xmin=-(n-1)*dx/2;

double x[n];
for (int i=0;i<n;i++) {x[i]=xmin+i*dx;}

double k[n];
for (int i=0;i<n;i++) {k[i]=2*M_PI/(n*dx)*(i-n/2);} 		

double data[2*n];
for (int i=0;i<n;i++)
{
if (x[i]==0) {data[2*i]=1;}
else {data[2*i]=sin(x[i])/x[i];}
data[2*i+1]=0;
}

gsl_fft_complex_radix2_forward(data,1,n);

double data_sort[2*n];
for (int i=0;i<n;i++)
{
if (i<n/2) {data_sort[2*i]=data[2*(i+n/2)]; data_sort[2*i+1]=data[2*(i+n/2)+1];}
if (i>=n/2) {data_sort[2*i]=data[2*(i-n/2)]; data_sort[2*i+1]=data[2*(i-n/2)+1];}
}

double f_t[n];
for (int i=0;i<n;i++)
{
double complex f_t_i=dx/sqrt(2*M_PI)*cexp(-I*k[i]*xmin)*(data_sort[2*i]+data_sort[2*i+1]*I);
f_t[i]=creal(f_t_i);
}

FILE *fp=NULL;
fp=fopen("ft_sinc2.txt","w");
for (int i=0;i<n;i++)
{fprintf(fp,"%e\t%e\n",k[i],f_t[i]);}

/* 
in linux terminal type (gnuplot must be installed):
gcc Q3_gsl_sinc.c -lgsl -lgslcblas -lm
./a.out
gnuplot
plot "ft_sinc2.txt" using 1:2 title "fourier transform of sinc(x) using gsl. ." with lines
set xlabel "k."
set ylabel "f(k)."
replot
*/

}
