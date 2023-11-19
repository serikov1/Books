#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>

//Решить задачу, записав уравнение Релея в виде системы ОДУ. Начальные
//условия:
double x_0 = 0;
double y_0 = 0.001;
double mu = 1000, Tk = 1000;
double h = 0.001;
double tau = -0.1;

double f(double x, double y) {
    return y;
}

double g (double x, double y) {
    return mu * (1 - std::pow(y, 2)) * y - x;
}

double K1(std::vector<double>& u) {
    double k1_0 = h * f(u[0], u[1]);
    double k1_n;
    while (true) {
        k1_n = k1_0 + tau * (h * (f(u[0] + k1_0 / std::sqrt(6), u[1]) - k1_0));
        if (std::abs(k1_n - k1_0) < 0.01 ) break;
        k1_0 = k1_n;
    }
    return  k1_n;
}

double L1(std::vector<double>& u){
    double l1_0 = h * g(u[0], u[1]);
    double l1_n;
    while (true) {
        l1_n = l1_0 + tau* (h * (g(u[0], u[1] + l1_0 / std::sqrt(6)) - l1_0));
        if (std::abs(l1_n - l1_0) < 0.01 ) break;
        l1_0 = l1_n;
    }
    return  l1_n;
}

double K2(std::vector<double>& u, double  k1, double l1){
    double k2_0 = h * f(u[0] + 1 * k1, u[1] + 1 * l1);
    double k2_n;
    while (true) {
        k2_n = k2_0 + tau* (h * (f(u[0] + 1 * k1 + 1 / std::sqrt(6) * k2_0, u[1] + 1 * l1) - k2_0));
        if (std::abs(k2_n - k2_0) < 0.01 ) break;
        k2_0 = k2_n;
    }
    return k2_n;
}

double L2(std::vector<double>& u, double  k1, double l1) {
    double l2_0 = h * g(u[0] + 1 * k1, u[1] + 1 * l1);
    double l2_n;
    while (true) {
        l2_n = l2_0 + tau* (h * (g(u[0] + 1 * k1, u[1] + 1 * l1 + 1 / std::sqrt(6) * l2_0) - l2_0));
        if (std::abs(l2_n - l2_0) < 0.001) break;
        l2_0 = l2_n;
    }
    return l2_n;
}
std::vector<double> delta(std::vector<double>& u) {
    double k1 = K1(u);
    double l1 = L1(u);
    double k2 = K2(u, k1, l1);
    double l2 = L2(u, k1, l1);
    return {(k1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * k2 ) / std::sqrt(6), (l1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * l2 ) / std::sqrt(6)};
}
std::vector<double> Runge_Kutta(std::vector<double>& u){
    u[0] += delta(u)[0];
    u[1] += delta(u)[1];
    return u;
}
std::vector<double> routine (std::vector<double>& U){
    uint32_t n = 0;
    std::vector<double> u(2);
    u = U;
    while(h * n <= Tk){
        u = Runge_Kutta(u);
        n++;
    }
    return u;
}
int main() {
    std::vector<double> u0 = {0, 0.001};
    std::vector<double> u = routine(u0);
    std::cout<< "x = " << u[0] << " y = " << u[1];
    return 0;
}
