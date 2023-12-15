#include <iostream>
#include "cmath"
#include "vector"
#include <cstdint>

double T = 100;
double h = 1;

double f (double x){
    return sin(x);
}

double K1(std::vector<double>& u){
    return  h * f(u[0]);
}

double K2(std::vector<double>& u){
    return h * f(u[0] + h / 10);
}

double K3(std::vector<double>& u){
    return h * f(u[0] + 2 * h / 10);
}

double K4(std::vector<double>& u){
    return h * f(u[0] + 3 * h / 10);
}

double K5(std::vector<double>& u){
    return h * f(u[0] + 8 * h / 10);
}

double K6(std::vector<double>& u){
    return h * f(u[0] + 9 * h / 10);
}

std::vector<double> delta(std::vector<double>& u) {
    double k1 = K1(u);
    double k2 = K2(u);
    double k3 = K3(u);
    double k4 = K4(u);
    double k5 = K5(u);
    double k6 = K6(u);
    return {h, (37 * k1 / 378 + 250 * k3 / 621 + 125 * k4 / 594 + 512 * k6 / 1771 )};
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