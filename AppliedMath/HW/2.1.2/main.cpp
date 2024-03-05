#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>
#include "linal.h"

double L = 10;
double gamma = 5./3;
double v_L = 0;
double rho_L = 13;
double P_L = 10;
double v_R = 0;
double rho_R = 1.3;
double P_R = 1;

double T = 0.02;
double tau = 10e-8;
uint32_t N = 1000;
double h = 10e-5;
int N_t = (int)(T / tau);
double lambda_max = 0;

double epsilon(double c)
{
    return pow(c, 2) / (gamma * (gamma - 1));
}

std::vector<std::vector<double>> omega_T(double u, double c, std::vector<std::vector<double>>& omega)
{
    omega[0][0] = -u * c;    omega[0][1] = c;     omega[0][2] = gamma - 1;
    omega[1][0] = -c * c;    omega[1][1] = 0;     omega[1][2] = gamma - 1;
    omega[2][0] = u * c;     omega[2][1] = -c;    omega[2][2] = gamma - 1;

    return omega;
}

std::vector<std::vector<double>> omega_T_rev(double u, double c, std::vector<std::vector<double>>& omega)
{
    omega[0][0] = 1 / (2 * pow(c, 2));            omega[0][1] = -2 / (2 * pow(c, 2));           omega[0][2] = 1 / (2 * pow(c, 2));
    omega[1][0] = (u + c) / (2 * pow(c, 2));      omega[1][1] = -2 * u / (2 * pow(c, 2));       omega[1][2] = (u - c) / (2 * pow(c, 2));
    omega[2][0] = 1 / (2 * (gamma - 1));                omega[2][1] = 0;                                    omega[2][2] = 1 / (2 * (gamma - 1));

    return omega;
}


std::vector<std::vector<double>>& lambda_mod(double u, double c, std::vector<std::vector<double>>& lambda)
{
    lambda[0][0] = std::abs(u + c);
    lambda[1][1] = std::abs(u);
    lambda[2][2] = std::abs(u - c);

    if (lambda[0][0] > lambda_max) lambda_max = lambda[0][0];
    if (lambda[1][1] > lambda_max) lambda_max = lambda[1][1];
    if (lambda[2][2] > lambda_max) lambda_max = lambda[2][2];

    return lambda;
}

std::vector<std::vector<double>>& A(double u, double c, std::vector<std::vector<double>>& a)
{
    a[0][0] = 0;                           a[0][1] = 1;                     a[0][2] = 0;
    a[1][0] = -pow(u, 2);            a[1][1] = 2 * u;                 a[1][2] = gamma - 1;
    a[2][0] = -gamma * u * epsilon(c);     a[2][1] = gamma * epsilon(c);    a[2][2] = u;

    return a;
}

std::vector<std::vector<double>> routine(std::vector<std::vector<double>>& w0, std::vector<std::vector<double>>& w){
    for (auto i = 1; i < N - 1; i++){
        double c = sqrt(gamma * (gamma - 1) * w0[i][2] / w0[i][0]);
        double u = w0[i][1] / w0[i][0];

        std::vector<std::vector<double>> A_i(3, std::vector<double>(3));
        A(u, c, A_i);

        std::vector<std::vector<double>> omega_T_i(3, std::vector<double>(3));
        omega_T(u, c, omega_T_i);

        std::vector<std::vector<double>> omega_T_rev_i(3, std::vector<double>(3));
        omega_T_rev(u, c, omega_T_rev_i);

        std::vector<std::vector<double>> lambda_mod_i(3, std::vector<double>(3, 0));
        lambda_mod(u, c, lambda_mod_i);

        std::vector<std::vector<double>> omega_lambda(3, std::vector<double>(3, 0));
        matrix_to_matrix(omega_T_rev_i, lambda_mod_i, omega_lambda);

        std::vector<std::vector<double>> Omega_Lambda_Omega(3, std::vector<double>(3, 0));
        matrix_to_matrix(omega_lambda, omega_T_i, Omega_Lambda_Omega);

        std::vector<double> two_w = {w0[i][0] * 2, w0[i][1] * 2, w0[i][2] * 2};

        std::vector<double> w_minus_two_w(3);
        vector_minus_vector(w0[i + 1], two_w, w_minus_two_w);

        std::vector<double> w_minus_w(3);
        vector_minus_vector(w0[i + 1], w0[i - 1], w_minus_w);

        std::vector<double> w_minus_two_w_plus_w(3);
        vector_plus_vector(w_minus_two_w, w0[i - 1], w_minus_two_w_plus_w);

        std::vector<double> A_i_to_w_minus_w(3);
        matrix_to_vector(A_i, w_minus_w, A_i_to_w_minus_w);

        std::vector<double> Omega_Lambda_Omega_w_minus_two_w_plus_w(3);
        matrix_to_vector(Omega_Lambda_Omega, w_minus_two_w_plus_w, Omega_Lambda_Omega_w_minus_two_w_plus_w);

        for (auto j = 0; j < 3; j++){
            w[i][j] = w0[i][j] - tau / (2 * h) * A_i_to_w_minus_w[j]
                      + tau / (2 * h) *  Omega_Lambda_Omega_w_minus_two_w_plus_w[j];
        }
    }
    w[0][0] = w[1][0];          w[0][1] = w[1][1];          w[0][2] = w[1][2];
    w[N - 1][0] = w[N - 2][0];  w[N - 1][1] = w[N - 2][1];  w[N - 1][2] = w[N - 2][2];

    return w;
}

std::vector<std::vector<std::vector<double>>> KER(std::vector<double>& v1, std::vector<double>& v2){
    std::vector<std::vector<std::vector<double>>> w0(N_t, std::vector<std::vector<double>>(N, std::vector<double>(3, 0)));
    std::vector<std::vector<std::vector<double>>> w(N_t, std::vector<std::vector<double>>(N, std::vector<double>(3, 0)));
    for(auto i = 0u; i < N / 2; i++) {
        w0[0][i][0] = v1[0];
        w0[0][i][1] = v1[0] * v1[1];
        w0[0][i][2] = v1[0] * v1[2];
    }
    for(auto i = N / 2u; i < N; i++) {
        w0[0][i][0] = v2[0];
        w0[0][i][1] = v2[0] * v2[1];
        w0[0][i][2] = v2[0] * v2[2];
    }
    routine(w0[0], w[0]);
    double cfl;
    for (auto i = 0u; i < N_t - 1; i++){
        routine(w[i], w[i + 1]);
        cfl = tau * lambda_max / h;
    }

    return w;
}

int main() {
    std::vector<double> v_l = {rho_L, v_L, P_L / ((gamma - 1) * rho_L)};
    std::vector<double> v_r = {rho_R, v_R, P_R / ((gamma - 1) * rho_R)};
    std::vector<std::vector<std::vector<double>>> solution = KER(v_l, v_r);

    std::cout<< "rho = " << std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << solution[150000][i][0] << ", " << std::endl;
    }
    std::cout<< "u = " << std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << solution[150000][i][1] / solution[150000][i][0] << ", " << std::endl;
    }
    std::cout<< "e = " << std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << solution[150000][i][2] / solution[150000][i][0] << ", " << std::endl;
    }
    std::cout<< "P = " << std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << solution[150000][i][2] * (gamma - 1) << ", " << std::endl;
    }

    return 0;
}
