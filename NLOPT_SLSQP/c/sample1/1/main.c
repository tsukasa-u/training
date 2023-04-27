#include <math.h>
#include <nlopt.h>
#include <stdio.h>

#include "weight.c"

int count = 0;
double myfunc(unsigned int n, const double *x, double *grad, void *my_func_data) {
    ++count;
    if (grad) {
        for (unsigned int j = 0; j < n - 1; ++j) grad[j] = 0.0;
        grad[n - 1] = 1.0;
    }
    printf("%f\n", x[n -1]);
    return x[n - 1];
}

typedef struct {
     double g;
     int k;
     unsigned int n;
     void* d;
} my_constraint_data;


double myconstraint_x(unsigned n, const double *x, double *grad, void *data) {
    my_constraint_data *d_p = (my_constraint_data *) data;
    double g = d_p->g;
    int k = d_p->k;
    int _n = d_p->n;
    double (*_d_p)[_n + 1][_n + 1] = (double (*)[_n + 1][_n + 1])d_p->d;
    if (grad) {
        unsigned int j = 0;
        for (; j < _n; ++j) grad[j] = (*_d_p)[k][j];
        for (; j < _n + _n + _n; ++j) grad[j] = 0.0;
        grad[n - 1] = 0.5*sqrt(2.0*g*x[_n + k])*cos(x[_n + _n + k]);
        grad[_n + k] = -(x[n - 1] - 0.0)/2.0*sqrt(2.0*g)*cos(x[_n + _n + k])/2.0/sqrt(x[_n + k]);
        grad[_n + _n + k] = (x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[_n + k])*sin(x[_n + _n + k]);
    }
    double sum = 0.0;
    for ( unsigned int j = 0; j <= _n; ++j) sum += (*_d_p)[k][j]*x[j];
    return sum - (x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[_n + k])*cos(x[_n + _n + k]);
}

double myconstraint_y(unsigned n, const double *x, double *grad, void *data) {
    my_constraint_data *d_p = (my_constraint_data *) data;
    double g = d_p->g;
    int k = d_p->k;
    int _n = d_p->n;
    double (*_d_p)[_n + 1][_n + 1] = (double (*)[_n + 1][_n + 1])d_p->d;
    if (grad) {
        unsigned int j = 0;
        for (; j < _n; ++j) grad[j] = 0.0;
        for (; j < _n + _n + _n; ++j) grad[j] = (*_d_p)[k][j - _n];
        grad[n - 1] = 0.5*sqrt(2.0*g*x[_n + k])*sin(x[_n + _n + k]);
        grad[_n + k] += -(x[n - 1] - 0.0)/2.0*sqrt(2.0*g)*sin(x[_n + _n + k])/2.0/sqrt(x[_n + k]);
        grad[_n + _n + k] = -(x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[_n + k])*cos(x[_n + _n + k]);
    }
    double sum = 0.0;
    for ( unsigned int j = 0; j <= _n; ++j) sum += (*_d_p)[k][j]*x[j + _n];
    return sum - (x[n - 1] - 0.0)/2.0*sqrt(2.0*g*x[_n + k])*sin(x[_n + _n + k]);
}

int main() {

    const unsigned int length_n = 10;
    const double _L = 0.5;

    // double w[10 + 1];
    // compute_weight_Chebyshev(w, 10, -0.5);

    double d[length_n + 1][length_n + 1];
    compute_D((void*)d, length_n, -0.5, -0.5);
        
    // double lb[3*(length_n + 1) + 1];
    // for (int i = 0; i < 3*(length_n + 1); ++i) lb[i]= -HUGE_VAL;
    // lb[3*(length_n + 1)] = 0;
    // nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
    nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, 3*(length_n + 1) + 1);
    
    
    // nlopt_set_lower_bounds(opt, lb);
    nlopt_set_lower_bounds1( opt, -HUGE_VAL);
    nlopt_set_upper_bounds1( opt, HUGE_VAL);

    nlopt_set_lower_bound( opt, 3*(length_n + 1), 0);
    nlopt_set_upper_bound( opt, 3*(length_n + 1), HUGE_VAL);

    {
        double lb[3*(length_n + 1) + 1], ub[3*(length_n + 1) + 1];
        nlopt_get_lower_bounds(opt, lb);
        nlopt_get_upper_bounds(opt, ub);

        for (int i = 0; i < 3*length_n + 4; ++i){
            printf("%f, %f\n", lb[i], ub[i]);
        }
    }

    nlopt_set_min_objective(opt, myfunc, NULL);

    my_constraint_data data[(length_n + 1)*2];
    for (unsigned int k = 0; k <= length_n; ++k) {
        data[k] = (my_constraint_data){1.0, k, length_n, (void*)d};
        nlopt_add_equality_constraint(opt, myconstraint_x, &data[k], 1e-8);
    }
    for (unsigned int k = 0; k <= length_n; ++k) {
        data[k + length_n] = (my_constraint_data){1.0, k, length_n, (void*)d};
        nlopt_add_equality_constraint(opt, myconstraint_y, &data[k + length_n + 1], 1e-8);
    }

    // nlopt_result nlopt_set_lower_bounds1(nlopt_opt opt, double lb);

    nlopt_set_xtol_rel(opt, 1e-4);
    double x[3*(length_n + 1) + 1];
    for (int i = 0; i <= 3*(length_n + 1); ++i) x[i] = 0.0;
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