#include <math.h>
#include <nlopt.h>
#include <stdio.h>

#include "weight.c"


typedef struct {
     double m_p, T;
     void* d;
     void* tau;
} my_constraint_data_x;

typedef struct {
     double u;
     unsigned int j;
} my_constraint_data_b;

int count = 0;
double myfunc(unsigned int n, const double *x, double *grad, void *f_data) {
    my_constraint_data_b *d_p = (my_constraint_data_b *) f_data;
    double u = d_p->u;
    unsigned int j = d_p->j;
    ++count;
    if (grad) {
        for (unsigned int k = 0; k < n; ++k) grad[k] = 0.0;
        grad[j] = -1.0;
    }
    for (unsigned int i = 0; i < n; ++i) printf("%f ", x[i]);
    printf("\n");
    return -x[j];
}

void myconstraint_r(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    double m_p = d_p->m_p;
    double T = d_p->T;
    double (*_d_p)[m][m] = (double (*)[m][m])d_p->d;
    double (*_tau_p)[m] = (double (*)[m])d_p->tau;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            for (unsigned int j =   0; j <   m; ++j) grad[k*n + j] = (*_d_p)[k][j];
            for (unsigned int j =   m; j < 5*m; ++j) grad[k*n + j] = 0.0;
            grad[k*n + 2*m + k] -= (x[n - 1] - x[n - 2])/2.0;
            grad[k*n + n - 1] = -0.5*x[2*m + k];
            grad[k*n + n - 2] =  0.5*x[2*m + k];
        }
    }
    for (unsigned int k = 0; k < m; ++k) {
        result[k] = 0.0;
        for ( unsigned int j = 0; j < m; ++j) result[j] += (*_d_p)[k][j]*x[j];
        result[2*m + k] += - (x[n - 1] - x[n - 2])/2.0*x[2*m + k];
    }

}

void myconstraint_theta(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    double m_p = d_p->m_p;
    double T = d_p->T;
    double (*_d_p)[m][m] = (double (*)[m][m])d_p->d;
    double (*_tau_p)[m] = (double (*)[m])d_p->tau;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            for (unsigned int j =   0; j <   m; ++j) grad[k*n + j] = 0.0;
            for (unsigned int j =   m; j < 2*m; ++j) grad[k*n + j] = (*_d_p)[k][j - m];
            for (unsigned int j = 2*m; j < 5*m; ++j) grad[k*n + j] = 0.0;
            grad[k*n       + k] -= -(x[n - 1] - x[n - 2])/2.0*x[3*m + k]/pow(x[k], 2.0);
            grad[k*n + 3*m + k] -= (x[n - 1] - x[n - 2])/2.0/x[k];
            grad[k*n + n - 1] = -0.5*x[3*m + k]/x[k];
            grad[k*n + n - 2] =  0.5*x[3*m + k]/x[k];
        }
    }
    for (unsigned int k = 0; k < m; ++k) {
        result[k] = 0.0;
        for ( unsigned int j = 0; j < m; ++j) result[j] += (*_d_p)[k][j]*x[m + j];
        result[2*m + k] += - (x[n - 1] - x[n - 2])/2.0*x[3*m + k]/x[k];
    }
}

void myconstraint_u(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    double m_p = d_p->m_p;
    double T = d_p->T;
    double (*_d_p)[m][m] = (double (*)[m][m])d_p->d;
    double (*_tau_p)[m] = (double (*)[m])d_p->tau;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            double t = (x[n - 1] - x[n - 2])/2.0*(*_tau_p)[k] + (x[n - 1] + x[n - 2])/2.0;
            for (unsigned int j =   0; j < 2*m; ++j) grad[k*n + j] = 0.0;
            for (unsigned int j = 2*m; j < 3*m; ++j) grad[k*n + j] = (*_d_p)[k][j - 2*m];
            for (unsigned int j = 3*m; j < 5*m; ++j) grad[k*n + j] = 0.0;
            grad[k*n       + k] -= (x[n - 1] - x[n - 2])/2.0*(-pow(x[3*m + k], 2.0)/pow(x[k], 2.0) + 2.0/pow(x[k], 2.0));
            grad[k*n + 3*m + k] -= (x[n - 1] - x[n - 2])/2.0*2.0*x[3*m + k]/x[k];
            grad[k*n + 4*m + k] -= (x[n - 1] - x[n - 2])/2.0*T*cos(x[4*m + k])/(1.0 - m_p*t);
            grad[k*n + n - 1] = -0.5*(pow(x[2*m + k], 2.0)/x[k] - 1.0/pow(x[k], 2.0) + T*sin(x[4*m + k])/(1.0 - m_p*t));
            grad[k*n + n - 2] =  0.5*(pow(x[2*m + k], 2.0)/x[k] - 1.0/pow(x[k], 2.0) + T*sin(x[4*m + k])/(1.0 - m_p*t));
        }
    }
    for (unsigned int k = 0; k < m; ++k) {
        double t = (x[n - 1] - x[n - 2])/2.0*(*_tau_p)[k] + (x[n - 1] + x[n - 2])/2.0;
        result[k] = 0.0;
        for ( unsigned int j = 0; j < m; ++j) result[j] += (*_d_p)[k][j]*x[m + j];
        result[2*m + k] += - (x[n - 1] - x[n - 2])/2.0*(pow(x[3*m + k], 2.0)/x[k] - 1.0/pow(x[k], 2.0) + T*sin(x[4*m + k])/(1.0 - m_p*t));
    }
}

void myconstraint_v(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_x *d_p = (my_constraint_data_x *) f_data;
    double m_p = d_p->m_p;
    double T = d_p->T;
    double (*_d_p)[m][m] = (double (*)[m][m])d_p->d;
    double (*_tau_p)[m] = (double (*)[m])d_p->tau;
    // grad[i*n + j]
    if (grad) {
        for (unsigned int k = 0; k < m; ++k) {
            double t = (x[n - 1] - x[n - 2])/2.0*(*_tau_p)[k] + (x[n - 1] + x[n - 2])/2.0;
            for (unsigned int j =   0; j < 2*m; ++j) grad[k*n + j] = 0.0;
            for (unsigned int j = 2*m; j < 3*m; ++j) grad[k*n + j] = (*_d_p)[k][j - 2*m];
            for (unsigned int j = 3*m; j < 5*m; ++j) grad[k*n + j] = 0.0;
            grad[k*n       + k] -= (x[n - 1] - x[n - 2])/2.0*x[2*m + k]*x[3*m + k]/pow(x[k], 2.0);
            grad[k*n + 2*m + k] -= -(x[n - 1] - x[n - 2])/2.0*x[3*m + k]/x[k];
            grad[k*n + 3*m + k] -= -(x[n - 1] - x[n - 2])/2.0*x[2*m + k]/x[k];
            grad[k*n + 4*m + k] -= -(x[n - 1] - x[n - 2])/2.0*T*sin(x[4*m + k])/(1.0 - m_p*t);
            grad[k*n + n - 1] = -0.5*(-x[3*m + k]*x[2*m + k]/x[k] + T*cos(x[4*m + k])/(1.0 - m_p*t));
            grad[k*n + n - 2] =  0.5*(-x[3*m + k]*x[2*m + k]/x[k] + T*cos(x[4*m + k])/(1.0 - m_p*t));
        }
    }
    for (unsigned int k = 0; k < m; ++k) {
        double t = (x[n - 1] - x[n - 2])/2.0*(*_tau_p)[k] + (x[n - 1] + x[n - 2])/2.0;
        result[k] = 0.0;
        for ( unsigned int j = 0; j < m; ++j) result[j] += (*_d_p)[k][j]*x[m + j];
        result[2*m + k] += - (x[n - 1] - x[n - 2])/2.0*(-x[3*m + k]*x[2*m + k]/x[k] + T*cos(x[4*m + k])/(1.0 - m_p*t));
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
    return x[j] - u;
}

double myconstraint_d(unsigned n, const double* x, double* grad, void* f_data) {
    my_constraint_data_b *d_p = (my_constraint_data_b *) f_data;
    double u = d_p->u;
    unsigned int j = d_p->j;
    if (grad) {
        for (unsigned int i = 0; i < n; ++i) grad[i] = 0.0;
        grad[j] = 2.0*x[j];
    }
    return pow(x[0], 2.0) - pow(u, 2.0);
}

int main() {

    const unsigned int length_n = 4;
    const double _L = 0.5;
    const unsigned int length_extra = 2;
    const unsigned int var_num = 5;

    double w[length_n + length_extra];
    compute_weight_Chebyshev(w, length_n);

    // for ( int i = 0; i <= length_n; ++i) printf("w%f\n", w[i]);
    // printf("\n");

    double d[length_n + 1][length_n + 1];
    double tau[length_n + 1];
    compute_D((void*)d, length_n, -0.5, -0.5, tau);
    
    for ( int i = 0; i <= length_n; ++i) {
        for ( int j = 0; j <= length_n; ++j) printf("%f ", d[i][j]);
        printf("\n");
    }

    // double lb[3*(length_n + 1) + 1];
    // for (int i = 0; i < 3*(length_n + 1); ++i) lb[i]= -HUGE_VAL;
    // lb[3*(length_n + 1)] = 0;
    // nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
    nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, var_num*(length_n + 1) + length_extra);
    
    
    // nlopt_set_lower_bounds(opt, lb);
    nlopt_set_lower_bounds1( opt, 0);
    nlopt_set_upper_bounds1( opt, HUGE_VAL);

    nlopt_set_lower_bound( opt, var_num*(length_n + 1), 0);
    nlopt_set_upper_bound( opt, var_num*(length_n + 1), HUGE_VAL);

    {
        double lb[var_num*(length_n + 1) + length_extra], ub[var_num*(length_n + 1) + length_extra];
        nlopt_get_lower_bounds(opt, lb);
        nlopt_get_upper_bounds(opt, ub);

        for (int i = 0; i < var_num*(length_n + 1) + length_extra; ++i){
            printf("%f, %f\n", lb[i], ub[i]);
        }
    }

    my_constraint_data_b data_myfunc = {length_n, 0};
    nlopt_set_min_objective(opt, myfunc, &data_myfunc);

    // double tol[var_num*(length_n + 1) + length_extra];
    // for (int i = 0; i < var_num*(length_n + 1) + length_extra; ++i) tol[i] = 1e-8;
    // my_constraint_data_x c_data_r = {1.0, 1.0, (void*) d, (void *)tau};
    // nlopt_add_equality_mconstraint( opt, length_n + 1, myconstraint_r, &c_data_r, tol);
    // my_constraint_data_x c_data_theta = {1.0, 1.0, (void*) d, (void *)tau};
    // nlopt_add_equality_mconstraint( opt, length_n + 1, myconstraint_theta, &c_data_r, tol);
    // my_constraint_data_x c_data_u = {1.0, 1.0, (void*) d, (void *)tau};
    // nlopt_add_equality_mconstraint( opt, length_n + 1, myconstraint_r, &c_data_u, tol);
    // my_constraint_data_x c_data_v = {1.0, 1.0, (void*) d, (void *)tau};
    // nlopt_add_equality_mconstraint( opt, length_n + 1, myconstraint_r, &c_data_v, tol);

    // nlopt_add_equality_constraint(opt, myconstraint_b, &(my_constraint_data_b){              0 , 1.0}, 1e-8);
    // nlopt_add_equality_constraint(opt, myconstraint_b, &(my_constraint_data_b){   length_n + 1 , 0.0}, 1e-8);
    // nlopt_add_equality_constraint(opt, myconstraint_b, &(my_constraint_data_b){2*(length_n + 1), 0.0}, 1e-8);
    // nlopt_add_equality_constraint(opt, myconstraint_b, &(my_constraint_data_b){3*(length_n + 1), 1.0}, 1e-8);

    my_constraint_data_b datas[length_n + 1];
    for (unsigned int i = (var_num - 1)*(length_n + 1); i < var_num*(length_n + 1); ++i) {
        datas[i - (var_num - 1)*(length_n + 1)] = (my_constraint_data_b){i, 1.0};
        nlopt_add_inequality_constraint(opt, myconstraint_d, &datas[i - (var_num - 1)*(length_n + 1)], 1e-4);
    }

    
    nlopt_set_ftol_abs( opt,  1e-4);
    nlopt_set_xtol_rel( opt,  1e-4);
    

    double x[var_num*(length_n + 1) + length_extra];
    for (int i = 0; i < var_num*(length_n + 1) + length_extra; ++i) x[i] = 1.0*((double)i/(double)length_n) + 1.0;
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