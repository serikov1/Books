#include <cstdint>
#include <iostream>
#include "cmath"
#include "vector"

double e = 2.718281828459045;
double l0 = e;
double T = pow(e, 2);
double hk = 0.001;
double hl = 0.001;
double y_l = 2 * pow(e, 2);


double h (double x, double y, double z) {
    return (1 / pow(x, 2) - exp(z) * y + e / log(x) * pow(y, 2));
}


double g (double x, double y, double z){
//    fmt::print("z = {}\n", z);
    return z;
}

double K1(std::vector<double>& u){
    return  hk * h(u[0], u[1], u[2]);
}

double L1(std::vector<double>& u){
    return  hl * g(u[0], u[1], u[2]);
}

double K2(std::vector<double>& u, double  k1, double l1){
    return hk * h(u[0] + 0.5 * hk, u[1] + 0.5 * l1, u[2] + 0.5 * k1);
}

double L2(std::vector<double>& u, double  k1, double l1) {
    return hl * g(u[0] + 0.5 * hl, u[1] + 0.5 * l1, u[2] + 0.5 * k1);
}

double K3(std::vector<double>& u, double  k2, double l2){
    return hk * h(u[0] + 0.5 * hk, u[1] + 0.5 * l2, u[2] + 0.5 * k2);
}

double L3(std::vector<double>& u, double  k2, double l2){
    return hl * g(u[0] + 0.5 * hl, u[1] + 0.5 * l2, u[2] + 0.5 * k2);
}

double K4(std::vector<double>& u, double  k3, double l3){
    return hk * h(u[0] + hk, u[1] + l3, u[2] + 0.5 * k3);
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

double t = 0.5;
std::vector<double> shoot(std::vector<double> U){
    uint16_t n = 0;
    double z = U[2];
    std::vector<double> u(3);
    u = U;

    while (std::abs(u[1] - y_l) > 0.001) {
        u = {e, e, z};
        while (l0 + hl * n < T) {
            u = calculate(u);
            n++;
        }

        z += t * (u[1] - y_l);
        n = 0;
    }
    //std::cout<<"z = " << z << std::endl;
    u = {e, e, z};
    while (l0 + hl * n < T) {
        u = calculate(u);
        if (n == 100) std::cout << "n = 100: x = "<< u[0] << " y = " << u[1] << std::endl;
        if (n == 500) std::cout << "n = 500: x = "<< u[0] << " y = " << u[1] << std::endl;
        if (n == 2000) std::cout << "n = 2000: x = "<< u[0] << " y = " << u[1] << std::endl;
        n++;
    }
    return u;
}

int main(){
    std::vector<double> u0 = {e, e, e};
    std::vector<double> u1 = shoot(u0);
    std::cout << "x = " << u1[0] << " y = " << u1[1] << std::endl;

    return 0;
}