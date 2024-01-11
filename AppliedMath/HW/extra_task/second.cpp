#include "cmath"
#include "vector"
#include <cstdint>
#include <iostream>

double T = 100;
double h = 0.1;

double f (double x){
    return sin(x);
}

double K1(std::vector<double>& u){
    return  h * f(u[0]);
}

double K2(std::vector<double>& u){
    return h * f(u[0] + h);
}

std::vector<double> delta(std::vector<double>& u) {
    double k1 = K1(u);
    double k2 = K2(u);
    return {h, (k1 / 2 + k2 / 2)};
}

std::vector<double> calculate(std::vector<double>& u){
    std::vector<double> d = delta(u);
    u[0] += d[0];
    u[1] += d[1];
    return u;
}

std::vector<double> routine (std::vector<double>& U){
    uint16_t n = 0;
    std::vector<double> u(2);
    u = U;
    while (h * n < T) {
        u = calculate(u);
        n++;
        std::cout << u[1] << ", " << std::endl;
    }
    return u;
}

int main(){
    std::vector<double> u_0 = {0, 1};
    std::vector<double> u = routine(u_0);
    return 0;
}