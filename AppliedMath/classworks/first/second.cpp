#include <iostream>
#include <cmath>
#include <vector>

double tau = 0.01;
double u_0 = 1;
double t_0 = 0;
double t_1 = 10;
double l = t_1/tau;
double u_i = u_0;
double k1(double u) {
    return u;
}
double k2(double k1, double x) {
    return x + k1*tau/2;
}
double k3(double k2, double x) {
    return x + k2 * tau/2;
}
double k4(double k3, double x) {
    return x + k3 * tau/2;
}
double u_i1;
int main() {
    for(auto i = 0; i < l; i++) {
        u_i1 = u_i + tau * 1/6 * (k1(u_i) + 2 * k2(u_i, k1(u_i)) + 2 * k3(u_i, k2(u_i, k1(u_i))) + k4(u_i, k3(u_i, k2(u_i, k1(u_i)))));
        u_i = u_i1;
        std::cout<<u_i1<<", " <<std::endl;
    }
}