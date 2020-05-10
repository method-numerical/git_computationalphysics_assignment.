#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include <complex.h> 

void main()
{
int n=40;
double dx=1.0;
double xmin=-(n-1)*dx/2;

double x[n];
double k[n];
fftw_complex fx[n],fk[n];
fftw_plan p;

for (int i=0;i<n;i++)
{
x[i]=xmin+i*dx;
k[i]=2*M_PI*(i-n/2)/(n*dx);
fx[i][0]=exp(-x[i]*x[i]);
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

double an_f_t[n];
for (int i=0;i<n;i++) {an_f_t[i]=exp(-k[i]*k[i]/4)/sqrt(2);}

FILE *f1=NULL;
f1=fopen("ft_gaussian.txt","w");
for (int i=0;i<n;i++)
{fprintf(f1,"%f\t%f\t%f\t%f\n",x[i],k[i],f_t[i],an_f_t[i]);} 

/*
then in  linux terminal type (gnuplot must be installed):
gcc Q4_fftw_gaussian.c -lfftw3 -lm
./a.out
gnuplot
plot "ft_gaussian.txt" using 2:3 title "fourier transform of gaussian . ." with lines lc rgbcolor "green","ft_gaussian.txt" u 2:4 title "analytical result . ." pt 1 ps 1
set xlabel "k."
set ylabel "f(k)."
replot
*/

}
