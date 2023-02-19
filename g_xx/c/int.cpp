#include <iostream>
#include <cstring>
#include "legendre.h"
#include <cmath>

//----------------

#define n 16
#define k 17

#define aH 0.0
#define bH 1.0

#define CH (bH - aH)/2.0

#define xH(t) ((bH - aH)*t/2.0 + (bH + aH)/2.0)

double (*P)(double) = P16;
double *x_k = x_17;

double WV[k];
double WH[k];

//----------------

const int k_max = 5;
double nabra_a[k_max][k_max];
double gamma_a[k_max];
double coeff_a[k_max][k_max];

double nabra_b[k_max][k_max];
double gamma_b[k_max];
double coeff_b[k_max][k_max];

const double epsiloon = 0.022;
const double L = 2.0;
const int N = 200;
double wall[3] = {0.0, 0.0, 0.0};
const double dx = L/N;
const double dt = 1E-7;
const int u_max = N;
const int steps_max = 10000000;

//----------------

// #define auto_eps 5E-8
#define auto_eps 1E0
#define auto_alpha 0.9

//----------------

void init() {
    for (int i = 0; i < k; ++i) {
        WH[i] = 2*(1 - x_k[i]*x_k[i])/(k*P(x_k[i]))/(k*P(x_k[i]));
    }
}

void gamma(double a[k_max], int l, int b) {

    auto f = [](double _n, int _k) {
        if (_k==0) return 1.0;
        double tmp = _n/_k;
        for (double i = 1.0; i < _k; ++i) {
            tmp *= (_n-i)/i;
        }
        return tmp;
    };

    for (int i = 0; i < l; ++i) {
        double sum = 0.0;
        for (int j = 0; j < k; ++j) {
            sum += WH[j]*f(b - xH(x_k[j]), i);
        }

        a[i] = (i%2==0) ? sum*CH : -1*sum*CH;
    }
    return;
}

void two_terms(double a[k_max][k_max], int l) {

    for (int i=0; i<l; ++i) for (int j=0; j<l; ++j) a[i][j] = 0.0;

    a[0][0] = 1.0;
    if (l==1) return;

    a[1][0] = a[1][1] = 1.0;
    if (l==2) return;

    for(int i = 2; i < l; ++i) {
        a[i][0] = 1.0;
        for (int j = 1; j<i; ++j) {
            a[i][j] = a[i-1][j-1] + a[i-1][j];
        }
        a[i][i] = 1.0;
    }
    return ;
}

void terms_k(double a[k_max][k_max], double b[k_max], int l) {

    for (int i=0; i<l; ++i) for (int j=0; j<l; ++j) a[i][j] = 0.0;

    a[0][0] = b[0]*nabra_a[0][0];
    if (l==1) return;

    for(int i = 1; i < l; ++i) {
        for (int j = 0; j <= i; ++j) {
            a[i][j] = a[i-1][j] + b[i]*nabra_a[i][j]*(j%2==0 ? 1: -1);
        }
    }
    return;
}

void update_wall(double u[u_max], double wall[4]) {
    wall[0] = u[u_max-2];
    wall[1] = u[u_max-1];
    wall[2] = u[0];
    wall[3] = u[1];
    return;
}

void init_u(double u[u_max]) {
    for (int i = 0; i < u_max; ++i) {
        u[i] = cos(2.0*M_PI*i/u_max);
    }
    return;
}

void cal_f(double u[u_max], double dx, double t, double f[u_max], double wall[3]) {
    int i = 0;
    f[  i] = (-u[i]*(u[i+1]-wall[1]) - epsiloon*epsiloon*(u[i+2]-2.0*u[i+1]+2.0*wall[1]-wall[0])/dx/dx)/2.0/dx;
    f[++i] = (-u[i]*(u[i+1]-u[i-1]) - epsiloon*epsiloon*(u[i+2]-2.0*u[i+1]+2.0*u[i-1]-wall[1])/dx/dx)/2.0/dx;
    for (i = 2;i < u_max-2; ++i) {
        f[i] = (-u[i]*(u[i+1]-u[i-1]) - epsiloon*epsiloon*(u[i+2]-2.0*u[i+1]+2.0*u[i-1]-u[i-2])/dx/dx)/2.0/dx;
    }
    f[  i] = (-u[i]*(u[i+1]-u[i-1]) - epsiloon*epsiloon*(wall[2]-2.0*u[i+1]+2.0*u[i-1]-u[i-2])/dx/dx)/2.0/dx;
    f[++i] = (-u[i]*(wall[2]-u[i-1]) - epsiloon*epsiloon*(wall[3]-2.0*wall[2]+2.0*u[i-1]-u[i-2])/dx/dx)/2.0/dx;
    return;
}

void liner(double u0[u_max], double u1[u_max], int l, double dt, double f[k_max][u_max], double b[k_max][k_max]) {
    double sum[u_max];
    for (int i = 0; i < u_max; ++i) sum[i] = 0.0;

    for (int i = 0; i < l; ++i) {
        // std::cout << ":" << b[l-1][i];
        for (int j = 0; j < u_max; ++j) {
            sum[j] += b[l - 1][i]*f[i][j];
        }
    }
    for (int i = 0; i < u_max; ++i) u1[i] = u0[i] + sum[i]*dt;
    return;
}

double calc_max(double u1[u_max], double u2[u_max]) {
    double tmp = 0.0;
    double max = 0.0;
    for (int i = 0; i < u_max; ++i) {
        tmp = std::abs(u2[2]-u1[1]);
        if (max < tmp) max = tmp;
    }
    return max;
}

void PECE(double u0[u_max], double u1[u_max], double u2[u_max], double *t, int l, double dt, double f0[k_max][u_max], double f1[k_max][u_max]) {

    
    liner(u0, u1, l, dt, f0, coeff_a);
    // std::cout << "t" << *t << " : dt" << dt << std::endl;
    
    update_wall(u1, wall);
    cal_f(u1, dx, *t, f1[0], wall);
    memcpy(&(f1[1]), f0, sizeof(f0[0])*(k_max-1));

    liner(u0, u2, l+1, dt, f1, coeff_b);

    update_wall(u2, wall);


    double beta = auto_alpha*std::pow((auto_eps/calc_max(u2, u1)), 1.0/l);
    // std::cout << "beta" << beta << std::endl;
    // if (beta > 1.0) {
        cal_f(u2, dx, *t, f1[0], wall);
    // } else {
    //     double ndiv = 10.0 - 1.0/beta;
    //     ndiv = ndiv <= 0.0 ? 10.0 : 10.0 - ndiv;

    //     for (int i = 0; i < ndiv; ++i) {
            // std::cout << "ndiv : " << ndiv << " : " << i << std::endl;
    //         PECE(u0, u1, u2, t, l, dt/ndiv, f0, f1);
    //         std::swap(u0, u2);
    //         std::swap(f0, f1);
    //     }
        
    // }

    *t += dt;
    return;
}

void output(int step, double t, double u[u_max]) {

    // std::cout << t;
    // for (int i = 0; i < u_max; ++i) std::cout << " " << u[i];
    // std::cout << std::endl;
    for (int i = 0; i < u_max; ++i) std::cout << t << "," << (double)L*i/N <<"," << u[i] << std::endl;

    return;
}

int main() {

    std::cout << "t,x,u" << std::endl;

    init();
    two_terms(nabra_a, k_max);

    gamma(gamma_a, k_max, 0);
    terms_k(coeff_a, gamma_a, k_max);

    gamma(gamma_b, k_max, 1);
    terms_k(coeff_b, gamma_b, k_max);

    // for (int i = 0; i < k_max; ++i) for (int j=0; j < k_max; ++j) std::cout << coeff_a[i][j] << " ";
    // std::cout << std::endl;
    // for (int i = 0; i < k_max; ++i) for (int j=0; j < k_max; ++j) std::cout << coeff_b[i][j] << " ";
    // std::cout << std::endl;

    double u0[u_max];
    double u1[u_max];
    double u2[u_max];
    double t = 0.0;

    init_u(u0);

    double f0[k_max][u_max];
    double f1[k_max][u_max];
    cal_f(u0, dx, 0.0, f0[0], wall);
    memset(f0[1], 0, sizeof(f0[0])*(k_max-1));

    // for (int i = 0; i < k_max; ++i)
    // memset(u0, 0, sizeof(u0));
    // memset(u1, 1, sizeof(u1));
    // for (int j = 0; j < u_max; ++j)
    // std::cout << u0[j] << " ";
    // std::cout << std::endl; 
    // for (int j = 0; j < u_max; ++j)
    // std::cout << u1[j] << " ";
    // std::cout << std::endl; 
    // memcpy(&(u1[1]), u0, sizeof(u0)-1);
    // for (int j = 0; j < u_max; ++j)
    // std::cout << u0[j] << " ";
    // std::cout << std::endl; 
    // for (int j = 0; j < u_max; ++j)
    // std::cout << u1[j] << " ";
    // std::cout << std::endl; 

    // u0[i] = i;
    // u0[++i] = i;
    // for (i = 2;i < u_max-2; ++i) {
    //     u0[i] = i;
    // }
    // u0[i] = i;
    // u0[++i] = i;

    // for (int j = 0; j < u_max; ++j)
    // std::cout << u0[j] << " ";
    // std::cout << std::endl; 


    for (int i = 0; i < steps_max; ++i) {
        PECE(u0, u1, u2, &t, 3, dt, f0, f1);
        std::swap(u0, u2);
        std::swap(f0, f1);
        if (i%10000==0) output(i, dt*i, u0);
    }

}