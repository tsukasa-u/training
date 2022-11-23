
#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

double f(double t, int N, double *a) {
    double tmp = 0.0;
    for (int i = 0; i <= N; ++i) {
        tmp += a[i]*std::pow(t, i);
    }
    return tmp;
}

double g(double t, int N, double *a) {
    double tmp = 0.0;
    for (int i = N; i >= 0; i--) {
        tmp = tmp*t + a[i];
    }
    return tmp;
}

void count(double (*func)(double, int, double*), double t, int N, double *a, double* result, long long* runtime) {

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    *result = (*func)(t, N, a);
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    *runtime = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
}

void ini_a(int N, double* a) {
    std::mt19937 mt(std::rand());
    std::uniform_int_distribution<> rand0_1(0.0, 1.0);
    for (int i = 0; i <= N; ++i) {
        a[i] = (double)rand0_1(mt);
    }
}

int main() {
    const int times[4] = {100000, 1000000, 10000000, 100000000};
    // const int times[4] = {1, 10, 100, 1000};
    double result_f[4] = {};
    long long runtime_f[4] = {};
    double result_g[4] = {};
    long long runtime_g[4] = {};

    const double t = 0.9;
    const int N_max = 100000000;
    double a_0;
    double* a = new double[N_max + 1];
    ini_a(N_max, a);

    for (int i = 0; i < 4; i++) {
        count(&f, t, times[i], a, &(result_f[i]), &(runtime_f[i]));
        count(&g, t, times[i], a, &(result_g[i]), &(runtime_g[i]));
    }

    for (int i = 0; i < 4; i++) std::cout << times[i] << " " << result_f[i] << " " << runtime_f[i] << std::endl;
    std::cout << std::endl;
    for (int i = 0; i < 4; i++) std::cout << times[i] << " " << result_g[i] << " " << runtime_g[i] << std::endl;
    
    free(a);
    return 1;
}