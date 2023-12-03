#include "cmath"
#include "vector"
#include <iostream>
#include <cstdint>

double e = 2.718281828459045;
double x0 = e;
double T = pow(e, 2);
double hk = 0.001;
double hl = 0.001;
double y_l = 2 * pow(e, 2);
double z0 = 2;

//    p =  0.5 * (-e^(log(x) + 1)) * x * log(x) / (sqrt(1 / x^2 + x^2 * e * log(x) - e^(log(x) + 1) * x * log(x))))
double p (double x){
    return 0.5 * (-pow(e, std::log(x) + 1)) * x * std::log(x) / (std::sqrt(1 / pow(x, 2) + pow(x, 2) * e * std::log(x) - pow(e, std::log(x) + 1) * x * log(x)));
}

//    q = x * e - 0.5 * e^(log(x) + 1)) / (sqrt(1 / x^2 + x^2 * e * log(x) - e^(log(x) + 1) * x * log(x))))
double q (double x){
    return (x * e - 0.5 * pow(e, log(x) + 1)) / (std::sqrt(1 / pow(x, 2) + pow(x, 2) * e * std::log(x) - pow(e, std::log(x) + 1) * x * log(x)));
}

//    r = sqrt(1 / x^2 + x^2 * e * log(x) - e^(log(x) + 1) * x * log(x))) - 1 / x
double r (double x){
    return (std::sqrt(1 / pow(x, 2) + pow(x, 2) * e * std::log(x) - pow(e, std::log(x) + 1) * x * log(x))) - 1 / x;
}

//    f =  p(x) * z + q(x) * y + r(x)
double f (double x, double y, double z){
    return p(x) * z + q(x) * y + r(x);
}

double h (double x, double y, double z) {
//    fmt::print("h = {} \n", 1 / pow(x, 2) - exp(z) * y + e / log(x) * pow(y, 2));
    return (1 / pow(x, 2) - exp(z) * y + e / log(x) * pow(y, 2));
}
double Y0 (double x) {
    return x * log(x);
}

//    z = g
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

std::vector<double> nu;
std::vector<double> Y;
std::vector<double> X;
std::vector<double> routine (std::vector<double>& U){
    uint16_t n = 0;
    std::vector<double> u(3);
    u = U;
    Y.push_back(u[1]);
    //std::cout<<"z = " << u[2] << std::endl;
    while (x0 + hl * n < T) {
        u = calculate(u);
        X.push_back(x0 + hl * n);
        nu.push_back(u[1]);
        Y.push_back(Y0(x0 + hl * n) + *(nu.end() - 1));

        if (n == 100) std::cout << "n = 100: x = "<< u[0] << " y = " << *(Y.end() - 1) << std::endl;
        if (n == 500) std::cout << "n = 500: x = "<< u[0] << " y = " << *(Y.end() - 1) << std::endl;
        if (n == 2000) std::cout << "n = 2000: x = "<< u[0] << " y = " << *(Y.end() - 1) << std::endl;
        n++;
    }

    return Y;
}


int main(){
    std::vector<double> u0 = {e, 0, z0};
    std::vector<double> u = routine(u0);
    std::cout << "x = " <<*(X.end() - 1)<< " y = "<< *(u.end() - 1)<< std::endl;
    return 0;
}