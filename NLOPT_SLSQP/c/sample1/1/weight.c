#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>

double compute_Chebyshev(unsigned int n, double x) {
    double r = sqrt(x*x - 1.0);
    return (pow(x + r, (double) n) + pow(x - r, (double) n))/2.0;
}

double compute_tau(unsigned int n, unsigned int k) {
    return cos((2.0*((double)k) + 1.0)*M_PI/(2.0*((double)k)));
}

void compute_weight_Chebyshev(double* w, unsigned int n, double alpha) {
    
// i = 0..(n+1)
    double double_alpha_p = pow(2.0, 2.0*alpha + 1.0);
    double N_2 = pow((double) n , 2.0);
    double gamma_N = gamma((double) n);
    double gamma_1_2 = pow(gamma((double) n + alpha + 1.0), 2.0);
    double gamma_2 = gamma((double) n + 2.0*alpha + 2.0);
    for (unsigned int i = 0; i <= n; ++i) {
        double x = compute_tau(n, i);
        double p_2 = pow(compute_Chebyshev(n, x), 2.0);
        w[i] = double_alpha_p*gamma_1_2/(N_2*gamma_N*gamma_2*p_2);
    }
}

void compute_D(void* d, unsigned int n, double alpha, double beta) {

    double (*d_p)[n + 1][n + 1] = (double (*)[n + 1][n + 1])d;

    double tau[n + 1];
    for (unsigned int k = 0; k <= n; ++k) tau[k] = compute_tau(n, k);

    double p_n[n + 1];
    for (unsigned int k = 0; k <= n; ++k) p_n[k] = compute_Chebyshev(n, tau[k]);

    (*d_p)[0][0] = (alpha - (double)n*((double)n + alpha + beta + 1.0))/2.0/(beta + 2.0);
    (*d_p)[n][n] = ((double)n*((double)n + alpha + beta + 1.0) - beta)/2.0/(alpha + 2.0);
    for (unsigned int k = 1; k < n - 1; ++k) for (unsigned int j = k + 1; j < n; ++j) (*d_p)[k][j] = p_n[k]/p_n[j]/(tau[k] - tau[j]);
    for (unsigned int k = 1; k < n; ++k) (*d_p)[k][k] = ((alpha + beta)*tau[k] + alpha - beta)/2.0/(1.0 - pow(tau[k], 2.0));
    for (unsigned int j = 1; j < n; ++j) (*d_p)[0][j] = -p_n[0]/p_n[j]/(beta + 1.0)/(1.0 + tau[j]);
    for (unsigned int j = 1; j < n; ++j) (*d_p)[n][j] = p_n[n]/p_n[j]/(alpha + 1.0)/(1.0 - tau[j]);
    for (unsigned int k = 1; k < n; ++k) (*d_p)[k][0] = p_n[k]/p_n[0]*(beta + 1.0)/(1 + tau[k]);
    for (unsigned int k = 1; k < n; ++k) (*d_p)[k][n] = -p_n[k]/p_n[n]*(alpha + 1.0)/(1 - tau[k]);
    (*d_p)[0][n] = -(alpha + 1.0)/2.0/(beta + 1.0)*p_n[0]/p_n[n];
    (*d_p)[n][0] = (beta + 1.0)/2.0/(alpha + 1.0)*p_n[n]/p_n[0];
}