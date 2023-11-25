#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>

// Решить краевую задачу
// y'' + (x^2 - 3)y' + (x^2-3)y*cos(x) = 2 - 6x + 2x^3 + (x^2 - 3)e^x*sin(x)(1 + cos(x)) + cos(x)(e^x + x^2 - 1 + x^4 - 3x^2)
// y(0) = 0
// y(pi) = pi^2

double epsilon = 0.0001;
double t = -0.005;
double x0 = 0;
double hk = 0.001;
double hl = 0.001;
double pi = 3.1415926535;
double T = pi;
double y_1 = pow(pi, 2);
double z0 = 1;

double f (double x, double y, double z) {
    return - (z + cos(x) * y) * (pow(x, 2) - 3) + 2 - 6 * x + 2 * pow(x, 3) + (pow(x, 2) - 3) * exp(x) + sin(x) *
    (1 + cos(x)) + cos(x) * (exp(x) + (pow(x, 2) - 1) + pow(x, 4) - 3 * pow(x, 2));
}

double g (double x, double y, double z){
    return z;
}

double K1(std::vector<double>& u){
    return  hk * f(u[0], u[1], u[2]);
}

double L1(std::vector<double>& u){
    return  hl * g(u[0], u[1], u[2]);
}

double K2(std::vector<double>& u, double  k1, double l1){
    return hk * f(u[0] + 0.5 * hk, u[1] + 0.5 * l1, u[2] + 0.5 * k1);
}

double L2(std::vector<double>& u, double  k1, double l1) {
    return hl * g(u[0] + 0.5 * hl, u[1] + 0.5 * l1, u[2] + 0.5 * k1);
}

double K3(std::vector<double>& u, double  k2, double l2){
    return hk * f(u[0] + 0.5 * hk, u[1] + 0.5 * l2, u[2] + 0.5 * k2);
}

double L3(std::vector<double>& u, double  k2, double l2){
    return hl * g(u[0] + 0.5 * hl, u[1] + 0.5 * l2, u[2] + 0.5 * k2);
}

double K4(std::vector<double>& u, double  k3, double l3){
    return hk * f(u[0] + hk, u[1] + l3, u[2] + 0.5 * k3);
}

double L4(std::vector<double>& u, double  k3, double l3){
    return hl * g(u[0] + hl, u[1] + l3, u[2] + 0.5 * k3);
}
std::vector<double> delta(std::vector<double>& u) {
    double k1 = K1(u);
    double l1 = L1(u);
    double k2 = K2(u, k1, l1);
    double l2 = L2(u, k1, l1);
    double k3 = K3(u, k2, l2);
    double l3 = L3(u, k2, l2);
    double k4 = K4(u, k3, l3);
    double l4 = L4(u, k3, l3);

    return {hl, (l1 + 2 * l2 + 2 * l3 + l4) / 6, (k1 + 2 * k2 + 2 * k3 + k4) / 6};
}
std::vector<double> calculate(std::vector<double>& u){
    std::vector<double> d = delta(u);
    u[0] += d[0];
    u[1] += d[1];
    u[2] += d[2];
    return u;
}

std::vector<double> getSolution(std::vector<double>& u) {
    std::vector<double> changeable(3);
    changeable = u;
    uint32_t n = 0;
    double z = u[2];
    while (std::abs(changeable[1] - y_1) > epsilon) {
        changeable = {0, 0, z};
        while (x0 + hl * n < T) {
            changeable = calculate(changeable);
            n++;
        }
        z += t * (changeable[1] - y_1);
        n = 0;
    }
    u[2] = z;
//        std::cout<<"z = " << z << std::endl;
    return changeable;
}

std::vector<double> routine (std::vector<double>& U){
    uint16_t n = 0;
    std::vector<double> u(3);
    u = U;
    while (x0 + hl * n < T) {
        u = calculate(u);
        if (hl * n == 0.5) std::cout<< "t = 0.5: " <<"x = " << u[0] << " y = " << u[1] << std::endl;
        if (hl * n == 1) std::cout<< "t = 1: " <<"x = " << u[0] << " y = " << u[1] << std::endl;
        if (hl * n == 1.5) std::cout<< "t = 1.5: " <<"x = " << u[0] << " y = " << u[1] << std::endl;
        if (hl * n == 2) std::cout<< "t = 2: " <<"x = " << u[0] << " y = " << u[1] << std::endl;
        if (hl * n == 2.5) std::cout<< "t = 2.5: " <<"x = " << u[0] << " y = " << u[1] << std::endl;
        n++;
    }
    return u;
}

int main() {
    std::vector<double> u0 = {0, 0, z0};
    std::vector<double> u = getSolution(u0);
    std::cout<<"x = " << u[0] << " y = " << u[1] << std::endl;
    u0 = {0, 0, u[2]};
    std::vector<double> u1 = routine(u0);
    return 0;
}
