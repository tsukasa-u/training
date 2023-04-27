#include <math.h>
#include <nlopt.h>
#include <stdio.h>

#include "weight.c"

int count = 0;
double myfunc(unsigned int n, const double *x, double *grad, void *my_func_data) {
    ++count;
    if (grad) {
        for (unsigned int j = 0; j < n; ++j) grad[j] = 1.0;
    }
    // for (int i = 0; i < n; ++i) printf("%f ", grad[i]);
    // printf("\n");
    // for (int i = 0; i < n; ++i) printf("%f ", x[i]);
    // printf("\n");
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += x[i];
    return sum;
}

typedef struct {
     double g;
     void* d;
} my_constraint_data_x;

typedef struct {
     double u;
     unsigned int j;
} my_constraint_data_b;

void myconstraint_x(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    // my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    // double g = d_p->g;
    // double (*_d_p)[m + 1][m + 1] = (double (*)[m + 1][m + 1])d_p->d;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            for (unsigned int j = 0; j < n; ++j) grad[k*n + j] = 0.0;
        }

        for (unsigned int k = 1; k < m - 1; ++k) {
            grad[k*n + k - 1] = cos(x[k - 1] - 2.0*x[k] + x[k + 1]);
            grad[k*n + k] = -2.0*cos(x[k - 1] - 2.0*x[k] + x[k + 1]);
            grad[k*n + k + 1] = cos(x[k - 1] - 2.0*x[k] + x[k + 1]);
        }
        {
            unsigned int k = 0;
            grad[k*n + n - 1] = cos(x[m - 1] - 2.0*x[k] + x[k + 1]);
            grad[k*n + k] = -2.0*cos(x[m - 1] - 2.0*x[k] + x[k + 1]);
            grad[k*n + k + 1] = cos(x[m - 1] - 2.0*x[k] + x[k + 1]);
            k = n - 1;
            grad[k*n + k - 1] = cos(x[k - 1] - 2.0*x[k] + x[0]);
            grad[k*n + k] = -2.0*cos(x[k - 1] - 2.0*x[k] + x[0]);
            grad[k*n] = cos(x[k - 1] - 2.0*x[k] + x[0]);
        }
    }
    {
        unsigned int k = 0;
        result[k*n + k] = cos(x[m - 1] - 2.0*x[k] + x[k + 1]);
        k = n - 1;
        result[k*n + k] = cos(x[k - 1] - 2.0*x[k] + x[0]);
    }
    for (unsigned int k = 1; k < m - 1; ++k) result[k] = sin(x[k - 1] - 2.0*x[k] + x[k + 1]);

}

void myconstraint_y(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    double g = d_p->g;
    double (*_d_p)[m + 1][m + 1] = (double (*)[m + 1][m + 1])d_p->d;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            for (unsigned int j = 0; j < m; ++j) grad[k*n + j] = 0.0;
            for (unsigned int j = m; j < 3*m; ++j) grad[k*n + j] = (*_d_p)[k][j-m];
            grad[k*n + n - 1] = 0.5*sqrt(2.0*g*x[m + k])*sin(x[m + m + k]);
            grad[k*n + m + k] += -(x[n - 1] - 0.0)/2.0*sqrt(2.0*g)*sin(x[m + m + k])/2.0/sqrt(x[m + k]);
            grad[k*n + m + m + k] = -(x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[m + k])*cos(x[m + m + k]);
        }
    }
    for (unsigned int k = 0; k < m; ++k) {
        result[k] = 0.0;
        for ( unsigned int j = 0; j <= m; ++j) result[j] += (*_d_p)[k][j]*x[j];
        result[k] += - (x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[m + k])*sin(x[m + m + k]);
    }

}

double myconstraint_b(unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_b *d_p = (my_constraint_data_b *) f_data;
    double u = d_p->u;
    unsigned int j = d_p->j;
    if (grad) {
        for (unsigned int i = 0; i < n; ++i) grad[i] = 0.0;
        grad[j] = 1.0;
    }
    return x[0] - u;
}

int main() {

    const unsigned int length_n = 800;
    const double _L = 0.5;

    // double w[length_n + 1];
    // compute_weight_Chebyshev(w, length_n);

    // for ( int i = 0; i <= length_n; ++i) printf("w%f\n", w[i]);
    // printf("\n");

    // double d[length_n + 1][length_n + 1];
    // compute_D((void*)d, length_n, -0.5, -0.5);
    
    // for ( int i = 0; i <= length_n; ++i) {
    //     for ( int j = 0; j <= length_n; ++j) printf("%f ", d[i][j]);
    //     printf("\n");
    // }

    // double lb[3*(length_n + 1) + 1];
    // for (int i = 0; i < 3*(length_n + 1); ++i) lb[i]= -HUGE_VAL;
    // lb[3*(length_n + 1)] = 0;
    // nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
    nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, length_n);
    
    
    // nlopt_set_lower_bounds(opt, lb);
    nlopt_set_lower_bounds1( opt, -1.0);
    nlopt_set_upper_bounds1( opt, HUGE_VAL);

    // nlopt_set_lower_bound( opt, length_n - 1, -10.0);
    // nlopt_set_upper_bound( opt, length_n - 1, HUGE_VAL);

    // {
    //     double lb[3*(length_n + 1) + 1], ub[3*(length_n + 1) + 1];
    //     nlopt_get_lower_bounds(opt, lb);
    //     nlopt_get_upper_bounds(opt, ub);

    //     for (int i = 0; i < 3*length_n + 4; ++i){
    //         printf("%f, %f\n", lb[i], ub[i]);
    //     }
    // }

    nlopt_set_min_objective(opt, myfunc, NULL);

    double tol[length_n];
    for (int i = 0; i < length_n; ++i) tol[i] = 1e-8;
    // my_constraint_data_x c_data = {1.0, (void*) d};
    nlopt_add_equality_mconstraint( opt, length_n, myconstraint_x, NULL, tol);

    // my_constraint_data_b data1 = {0, 0.0};
    // nlopt_add_equality_constraint(opt, myconstraint_b, &data1, 1e-8);

    // my_constraint_data_b data2 = {length_n + 1, 0.0};
    // nlopt_add_equality_constraint(opt, myconstraint_b, &data2, 1e-8);

    // my_constraint_data_b data3 = {length_n, _L};
    // nlopt_add_equality_constraint(opt, myconstraint_b, &data3, 1e-8);

    nlopt_set_ftol_abs( opt,  1e-4);
    nlopt_set_xtol_rel( opt,  1e-4);

    double x[length_n];
    for (int i = 0; i < length_n; ++i) x[i] = 1.0*i;
    double minf;
    nlopt_result result = nlopt_optimize(opt, x, &minf);
    if (result < 0) {
        printf("nlopt failed!\n");
        printf("err: %d\n", result);
    }
    else {
        printf("found minimum after %d evaluations\n", count);
        printf("found minimum at f = %0.10g\n", minf);
    }

    nlopt_destroy(opt);
    return 1;
}