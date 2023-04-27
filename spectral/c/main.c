#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* #include "sslib.h" */
/* #include "pssub.h" */
#include "fft.c"

#define N  256
#define ITEL 8
// #define NT 31847*16
// #define NTINT 31847

#define NT 6000*128
#define NTINT 6000
// #define N  4
// #define ITEL 2
// #define NT 1*2
// #define NTINT 1

void initc(double [],double [],double);
void advance(double [],double [],double,double,double);
void rk(double [],double [],double,double,double);
void rhs(double [],double [],double [],double [],double,double);
void fft1(double *, double *, int, int, int);

int main()
{
     double u[N],ar[N],ai[N];
     double vr[N],vi[N];
     double dx,dt,delta;
	 int i,j;

     dx=2.0/(double)(N);
//      dt=0.00001;
dt = 0.1;
     delta=0.022*0.022;

     printf("t,x,u\n"); 
            

     i=0;
     initc(ar,ai,dx);

            for(j=0;j<N;j++) {vr[j]=ar[j];vi[j]=ai[j];}
            fft1(vr,vi,N,ITEL,1);
            for(j=0;j<N;j++)  printf("%d, %f, %f\n", i, (double)j*dx,vr[j]); 
            printf("\n\n");

     for (i=1;i<=NT;i++) {
        advance(ar,ai,delta,dx,dt);
        if((i%NTINT)==0) {

            for(j=0;j<N;j++) {vr[j]=ar[j];vi[j]=ai[j];}
            fft1(vr,vi,N,ITEL,1);
            for(j=0;j<N;j++)  printf("%d, %f, %f\n", i, (double)j*dx,vr[j]); 
            printf("\n\n");

        }
     }

     return 0;
}

void initc(double ar[],double ai[],double dx)
{
/*
      ar, ai in wave space
*/
      int i;double u[N],pi;
      pi=4.0*atan(1.0);
      for (i=0;i<N;i++) {u[i]=cos(pi*dx*i); /* printf("%f %f\n",dx*i,u[i]); */}
      for (i=0;i<N;i++) {ar[i]=u[i];ai[i]=0.0;}
      fft1(ar,ai,N,ITEL,0);
      return;
}

void  advance(double ar[],double ai[],double delta,double dx,double dt)
{
      rk(ar,ai,delta,dx,dt);
      return;
}

void rk(double ar[], double ai[], double delta, double dx, double dt)
{
     int i;
     double aro[N],aio[N];
     double ak1r[N],ak1i[N],ak2r[N],ak2i[N],ak3r[N],ak3i[N],ak4r[N],ak4i[N];
     for (i=0;i<N;i++)  {aro[i]=ar[i];aio[i]=ai[i];}
     rhs(ar,ai,ak1r,ak1i,delta,dt);
     for (i=0;i<N;i++)  {ar[i]=aro[i]+ak1r[i]*0.5;ai[i]=aio[i]+ak1i[i]*0.5;}
     rhs(ar,ai,ak2r,ak2i,delta,dt);
     for (i=0;i<N;i++)  {ar[i]=aro[i]+ak2r[i]*0.5;ai[i]=aio[i]+ak2i[i]*0.5;}
     rhs(ar,ai,ak3r,ak3i,delta,dt);
     for (i=0;i<N;i++)  {ar[i]=aro[i]+ak3r[i];ai[i]=aio[i]+ak3i[i];}
     rhs(ar,ai,ak4r,ak4i,delta,dt);
     for (i=0;i<N;i++)  {ar[i]=aro[i]+(ak1r[i]+2.0*ak2r[i]+2.0*ak3r[i]+ak4r[i])/6.0;
                         ai[i]=aio[i]+(ak1i[i]+2.0*ak2i[i]+2.0*ak3i[i]+ak4i[i])/6.0;}
     return;
}

void rhs(double ar[],double ai[],double akr[],double aki[],double delta,double dt)
{
/* INPUT  : ar , ai  in wave number space  */
/* OUTPUT : akr aki  = rhs *dt in wave number space  */
     int i;
     double ur[N],ui[N],dur[N],dui[N];
     double pi,di,d3i;

     double dx;
     dx=2.0/(double)(N);

     pi=atan(1.0)*4.0;
     for(i=0;i<N;i++) {ur[i]=ar[i];ui[i]=ai[i];}
     fft1(ur,ui,N,ITEL,1);

     for(i=0;i<=N/2;i++) {di=(double)i;
	                dur[i]=-1.0*pi*di*ai[i];
                        dui[i]= 1.0*pi*di*ar[i];
                         }
     for(i=1;i<N/2;i++) {dur[N-i]=dur[i];
                         dui[N-i]=-dui[i];}
     fft1(dur,dui,N,ITEL,1);

     for(i=0;i<N;i++) {dur[i]=ur[i]*dur[i];dui[i]=ui[i]*dui[i];}

     fft1(dur,dui,N,ITEL,0);
     for (i=0;i<=N/2;i++) {di=(double)i;
                     d3i=di*di*di*pi*pi*pi;
                     akr[i]=-dt*(dur[i]+1.0*d3i*delta*ai[i]);
                     aki[i]=-dt*(dui[i]-1.0*d3i*delta*ar[i]);
                      }
     for (i=1;i<N/2;i++) {di=(double)i;
                     akr[N-i]=akr[i];
                     aki[N-i]=-aki[i];
                      }

     return;
}
