#include "cmath"
#include "vector"
#include <iostream>
#include <fstream>
#include <cstdint>

//uint8_t B = 50;
uint16_t L = 500;
uint8_t p_start = 100;
uint8_t p_inj = 150;
uint8_t p_prod = 50;
double k = 10e-12;
double mu = 10e-9;
double fi = 0.2;
double c_f = 10e-5;
double rho0 = 1000;
double p0 = 120;
uint16_t h = 1;
uint16_t tau = 1;
uint16_t T = 240;

double rho(double p){
    return rho0 * (1 + c_f * (p - p0));
}

double C(double rho){
    return (k * rho) / (mu * h);
}

double B(double rho) {
    return (k * rho) / (mu * h);
}

double A(double c, double b){
    return - c - b - fi * c_f * rho0 / tau;
}

double D(double p){
    return - fi * c_f * rho0 * p / tau;
}

void run_through(std::vector<double> &p, std::vector<double> &p1){
    uint16_t N = p.size();

    std::vector<double> a(N);
    std::vector<double> b(N);
    std::vector<double> c(N);
    std::vector<double> d(N);
    std::vector<double> alpha(N);
    std::vector<double> betta(N);

    a[0] = 1;
    b[0] = 0;
    a[N - 1] = 1;
    c[N - 1] = 0;
    d[0] = p_inj;
    d[N - 1] = p_prod;
    alpha[0] = - b[0] / a[0];
    betta[0] = d[0] / a[0];

    for (auto i = 1u; i < N - 1; i++){
        double p_plus = p[i] >= p[i + 1] ? p[i] : p[i + 1];
        double p_minus = p[i - 1] >= p[i] ? p[i - 1] : p[i];
        c[i] = C(rho(p_minus));
        b[i] = B(rho(p_plus));
        a[i] = A(c[i], b[i]);
        d[i] = D(p[i]);

        alpha[i] = - b[i - 1] / (c[i - 1] * alpha[i - 1] + a[i - 1]);
        betta[i] = (d[i - 1] - c[i - 1] * betta[i - 1]) / (c[i - 1] * alpha[i - 1] + a[i - 1]);
    }
    alpha[N - 1] = - b[N - 2] / (c[N - 2] * alpha[N - 2] + a[N - 2]);
    betta[N - 1] = (d[N - 2] - c[N - 2] * betta[N - 2]) / (c[N - 2] * alpha[N - 2] + a[N - 2]);

    p1[N - 1] = (d[N - 1] - c[N - 1] * betta[N - 1]) / (c[N - 1] * alpha[N - 1] + a[N - 1]);
    for (auto i = N - 2; i >= 0; i--) {
        p1[i] = alpha[i + 1] * p1[i + 1] + betta[i + 1];
    }
}

std::vector<double> routine(std::vector<double>& p){
    std::vector<double> p1(p.size());
    std::vector<double> p2(p.size());
    std::ofstream out;
    out.open("preasure.txt");
    run_through(p, p1);
    for (auto i = 1u; i < T / tau; i++){
        run_through(p1, p2);
        std::copy(p2.begin(), p2.end(), p1.begin());

        if (out.is_open()) {
//            out <<  "\n t = " << i << ":" << std::endl;
            for (auto j = 0u; j < L / h ; j++){
                out <<  p1[j] << ", ";
            }
        }
    }
    out.close();
    return p1;
}

int main(){
    std::vector<double> p(L / h, p_start);
    p[0] = p_inj;
    p[L / h - 1] = p_prod;

    std::vector<double> p1 = routine(p);
    for (auto i = 0u; i < L / h ; i++){
       std::cout<<"p_i = " << p1[i] << std::endl;
    }

    return 0;
}
