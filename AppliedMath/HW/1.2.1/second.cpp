#include "cmath"
#include "vector"
#include <fmt/core.h>

double T = 1;
double hk = 0.001;
double hl = 0.001;
double y_l = 2;


// dz = f(x, y, z) * dx, where f(x, y, z) = x * sqrt(y)
double f (double x, double y, double z){
    return x * std::sqrt(y);
}

// dy = g(x, y, z) * dx, where g(x, y, z) = z
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
//    fmt::print("1 = {}, 2 = {}\n", Theta1(k1, k2, k3), Theta2(l1, l2, l3));
    return {hl, (l1 + 2 * l2 + 2 * l3 + l4) / 6, (k1 + 2 * k2 + 2 * k3 + k4) / 6};
}



std::vector<double> Runge_Kutte(std::vector<double>& u){
    std::vector<double> d = delta(u);
    u[0] += d[0];
    u[1] += d[1];
    u[2] += d[2];
    return u;
}


std::vector<double> routine (std::vector<double>& U){
    uint16_t n = 0;
    uint16_t k = 1;
    double z = U[2];
    std::vector<double> u(3);
    u = U;

    while (std::abs(u[1] - y_l) > 0.001) {
        u = {0, 0, z};
        while (hl * n < T) {
            u = Runge_Kutte(u);
            n++;
        }
//        fmt::print("y = {}\n", u[1]);
        if (u[1] > y_l) z -= std::abs(u[1] - y_l) ;
        else z += std::abs(u[1] - y_l);
        n = 0;
        k++;
    }
    fmt::print("z = {}\n", z);
    return u;
}


int main(){
    double z0 = 2;
    std::vector<double> u0 = {0, 0, z0};
    std::vector<double> u = routine(u0);
    fmt::print("x = {}, y = {}\n",u[0], u[1]);
    return 0;
}
