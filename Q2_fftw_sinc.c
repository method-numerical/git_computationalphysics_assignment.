#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include <complex.h> 

void main () 
{
int n=512;
double dx=1.0;
double xmin=-(n-1)*dx/2;
double x[n], k[n];

fftw_complex fx[n], fk[n];
fftw_plan p;

for (int i=0;i<n;i++)
{
x[i]=xmin+i*dx;
k[i]=2*M_PI/(n*dx)*(i-n/2);
if (x[i]==0) {fx[i][0]=1;} else {fx[i][0]=sin(x[i])/x[i];}
fx[i][1]=0;
}

p=fftw_plan_dft_1d(n,fx,fk,FFTW_FORWARD,FFTW_ESTIMATE);
fftw_execute(p); 

fftw_complex fk_sort[n];
for (int i=0;i<n;i++)
{
if (i<n/2) {fk_sort[i][0]=fk[i+n/2][0]; fk_sort[i][1]=fk[i+n/2][1];}
if (i>=n/2) {fk_sort[i][0]=fk[i-n/2][0]; fk_sort[i][1]=fk[i-n/2][1];}
}

double f_t[n];
for (int i=0;i<n;i++)
{
double complex f_t_i=dx*sqrt(n/(2*M_PI))*cexp(-I*k[i]*xmin)*(fk_sort[i][0]+fk_sort[i][1]*I)/sqrt(n);
f_t[i]=creal(f_t_i);
}

FILE *fp=NULL;
fp=fopen("ft_sinc1.txt","w");
for (int i=0;i<n;i++)
{fprintf(fp,"%e\t%e\n",k[i],f_t[i]);}

/*
then in  linux terminal type (gnuplot must be installed): 
gcc Q2_fftw_sinc.c -lfftw3 -lm
./a.out
gnuplot
plot "ft_sinc1.txt" using 1:2 title "fourier transform of sinc(x) using fftw . ." with lines
set xlabel "k."
set ylabel "f(k)."
replot
*/

}
