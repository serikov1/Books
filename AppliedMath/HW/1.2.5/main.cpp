#include <iostream>
#include <cstdint>
#include "cmath"
#include "vector"

double T = 17.0652165601579625588917206249;
double T0 = 0;
double h = 0.001;
double nu = 0.012277471;
double eta = 1 - nu;
double tau = -0.0000001;

// dx = f1(x, y) * dt
double f1 (double x, double y, double u, double z){
    return u;
}

// dy = f2(x, y) * dt
double f2 (double x, double y, double u, double z){
    return z;
}

// du = f3(x, y) * dt
double f3 (double x, double y, double u, double z){
    double A = std::sqrt(std::pow(std::pow(x + nu, 2) + std::pow(y, 2), 3));
    double B = std::sqrt(std::pow(std::pow(x - eta, 2) + std::pow(y, 2), 3));
    return x + 2 * z - eta * (x + nu) / A - nu * (x - eta) / B;
}


// dz = f4(x, y) * dt
double f4 (double x, double y, double u, double z){
    double A = std::sqrt(std::pow(std::pow(x + nu, 2) + std::pow(y, 2), 3));
    double B = std::sqrt(std::pow(std::pow(x - eta, 2) + std::pow(y, 2), 3));
    return y - 2 * u - eta * y / A - nu * y / B;
}


double K1(std::vector<double>& u){
    double k1_0 = 0;
    double k1_n;
    while (true) {
        k1_n = k1_0 + tau * (h * (f1(u[0] + k1_0 / std::sqrt(6), u[1], u[2], u[3]) - k1_0));
        if (std::abs(k1_n - k1_0) < 0.01 ) break;
        k1_0 = k1_n;
    }
    return  k1_n;
}

double L1(std::vector<double>& u){
    double l1_0 = 0;
    double l1_n;
    while (true) {
        l1_n = l1_0 + tau * (h * (f2(u[0], u[1] + l1_0 / std::sqrt(6), u[2], u[3]) - l1_0));
        if (std::abs(l1_n - l1_0) < 0.01 ) break;
        l1_0 = l1_n;
    }
    return  l1_n;
}

double P1(std::vector<double>& u){
    double p1_0 = 0;
    double p1_n;
    while (true) {
        p1_n = p1_0 + tau * (h * (f3(u[0], u[1], u[2] + p1_0 / std::sqrt(6), u[3]) - p1_0));
        if (std::abs(p1_n - p1_0) < 0.01 ) break;
        p1_0 = p1_n;
    }
    return  p1_n;
}

double M1(std::vector<double>& u){
    double m1_0 = 0;
    double m1_n;
    while (true) {
        m1_n = m1_0 + tau * (h * (f4(u[0], u[1], u[2], u[3] + m1_0 / std::sqrt(6)) - m1_0));
        if (std::abs(m1_n - m1_0) < 0.01 ) break;
        m1_0 = m1_n;
    }
    return  m1_n;
}

double K2(std::vector<double>& u, double  k1, double l1, double p1, double m1){
    double k2_0 = 0;
    double k2_n;
    while (true) {
        k2_n = k2_0 + tau * (h * (f1(u[0] + 1 * k1 + 1 / std::sqrt(6) * k2_0, u[1] + 1 * l1, u[2] + 1 * p1, u[3] + 1 * m1) - k2_0));
        if (std::abs(k2_n - k2_0) < 0.001 ) break;
        k2_0 = k2_n;
    }
    return k2_n;
}

double L2(std::vector<double>& u, double  k1, double l1, double p1, double m1) {
    double l2_0 = 0;
    double l2_n;
    while (true) {
        l2_n = l2_0 + tau * (h * (f2(u[0] + 1 * k1, u[1] + 1 * l1 + 1 / std::sqrt(6) * l2_0, u[2] + 1 * p1, u[3] + 1 * m1) - l2_0));
        if (std::abs(l2_n - l2_0) < 0.001) break;
        l2_0 = l2_n;
    }
    return l2_n;
}

double P2(std::vector<double>& u, double  k1, double l1, double p1, double m1){
    double p2_0 = 0;
    double p2_n;
    while (true) {
        p2_n = p2_0 + tau * (h * (f3(u[0] + 1 * k1, u[1] + 1 * l1, u[2] + 1 * p1 + 1 / std::sqrt(6) * p2_0, u[3] + 1 * m1) - p2_0));
        if (std::abs(p2_n - p2_0) < 0.001 ) break;
        p2_0 = p2_n;
    }
    return p2_n;
}

double M2(std::vector<double>& u, double  k1, double l1, double p1, double m1) {
    double m2_0 = 0;
    double m2_n;
    while (true) {
        m2_n = m2_0 + tau * (h * (f2(u[0] + 1 * k1, u[1] + 1 * l1, u[2] + 1 * p1, u[3] + 1 * m1 + 1 / std::sqrt(6) * m2_0) - m2_0));
        if (std::abs(m2_n - m2_0) < 0.001) break;
        m2_0 = m2_n;
    }
    return m2_n;
}

std::vector<double> delta(std::vector<double>& u) {
    double k1 = K1(u);
    double l1 = L1(u);
    double p1 = P1(u);
    double m1 = M1(u);
    double k2 = K2(u, k1, l1, p1, m1);
    double l2 = L2(u, k1, l1, p1, m1);
    double p2 = P2(u, k1, l1, p1, m1);
    double m2 = M2(u, k1, l1, p1, m1);
    return {(k1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * k2 ) / std::sqrt(6), (l1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * l2 ) / std::sqrt(6),
            (p1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * p2 ) / std::sqrt(6), (m1 * (3 + std::sqrt(6)) + (3 - std::sqrt(6)) * m2 ) / std::sqrt(6)};
}



std::vector<double> calculate(std::vector<double>& u){
    std::vector<double> d = delta(u);
    u[0] += d[0];
    u[1] += d[1];
    u[2] += d[2];
    u[3] += d[3];
    return u;
}



std::vector<double> routine (std::vector<double>& U){
    uint32_t n = 0;
    std::vector<double> u(4);
    u = U;
    while(T0 + h * n <= T){
        u = calculate(u);
        n++;
    }
    return u;
}

int main(){
    std::vector<double> u0 = {0.994, 0, 0, -2.00158510637908252240537862224}; // начальные условия
    std::vector<double> u = routine(u0);
    std::cout<< "T = 1 : x = "<< u[0] << " y = " << u[1] << " u = " << u[2] << " z = " << u[3]<< std::endl;

    T *= 50;
    std::vector<double> u2 = routine(u0);
    std::cout<< "T = 50 : x = "<< u[0] << " y = " << u[1] << " u = " << u[2] << " z = " << u[3]<< std::endl;

    T *= 2;
    std::vector<double> u1 = routine(u0);
    std::cout<< "T = 100 : x = "<< u[0] << " y = " << u[1] << " u = " << u[2] << " z = " << u[3]<< std::endl;
    return 0;
}