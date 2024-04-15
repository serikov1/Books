#include "cmath"
#include "vector"
#include <fmt/core.h>
#include <fstream>

#define PI 3.1415

uint8_t Lx = 1;
uint8_t Ly = 1;
uint8_t T = 1;
double lambda = 10e-5;

double hx = 0.01;
double hy = 0.01;
double tau = 0.01;
double rx = 25 * lambda * tau / pow(hx, 2);
double ry = lambda * tau / pow(hy, 2);

uint16_t Nx = uint16_t (double(Lx) / hx);
uint16_t Ny = uint16_t (double(Ly) / hy);

double f_t0(double x, double y){
    return cos(PI * x) * sin(5 * PI * y);
}

double f_x0(double t, double y){
    return sin(5 * PI * y) * exp(-50 * pow(PI, 2) * lambda * t);
}

double f_x1(double t, double y){
    return -sin(5 * PI * y) * exp(-50 * pow(PI, 2) * lambda * t);
}

double f0(double t, double x, double y){
    return cos(PI * x) * sin(5 * PI * y) * exp(- 50 * pow(PI, 2) * lambda * t);
}

void run_through_x(std::vector<std::vector<double>> &f, std::vector<std::vector<double>> &f1, double t){
    std::vector<double> a(Nx);
    std::vector<double> b(Nx);
    std::vector<double> c(Nx);
    std::vector<std::vector<double>> d(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> alpha(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> betta(Nx, std::vector<double>(Ny));

    a[0] = 1;
    b[0] = 0;
    c[0] = 0;
    for (auto i = 1u; i < Nx - 1; i++){
        a[i] = 2 * (1 + rx) / tau;
        b[i] = rx / tau;
        c[i] = rx / tau;
    }
    a[Nx - 1] = 1;
    b[Nx - 1] = 0;
    c[Nx - 1] = 0;

    for (auto j = 0u; j < Ny; j++){
        d[0][j] = f_x0(t, hy * j);
        d[Nx - 1][j] = f_x1(t, hy * j);
    }

    for (auto j = 0u; j < Ny; j++){
        for (auto i = 1u; i < Nx - 1; i++) {
            d[i][j] = 2 * f[i][j] / tau;
        }
    }

    for (auto j = 0u; j < Ny; j++){
        alpha[0][j] = b[0] / a[0];
        betta[0][j] = d[0][j] / a[0];
        for (auto i = 1u; i < Nx; i++) {
            alpha[i][j] = b[i] / (- c[i] * alpha[i - 1][j] + a[i]);
            betta[i][j] = (d[i][j] + c[i] * betta[i - 1][j]) / (- c[i] * alpha[i - 1][j] + a[i]);
        }
    }
    for (auto j = Ny - 1; j >= 0; j--){
        f1[Nx - 1][j] = betta[Nx - 1][j];
        for (auto i = Nx - 2; i >= 0; i--){
            f1[i][j] = alpha[i][j] * f1[i + 1][j] + betta[i][j];
        }
    }

}

void run_through_y(std::vector<std::vector<double>> &f, std::vector<std::vector<double>> &f1, double t){

    std::vector<double> a(Ny);
    std::vector<double> b(Ny);
    std::vector<double> c(Ny);
    std::vector<std::vector<double>> d(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> alpha(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> betta(Nx, std::vector<double>(Ny));

    a[0] = 1;
    b[0] = 0;
    c[0] = 0;
    for (auto i = 1u; i < Nx - 1; i++){
        a[i] = 2 * (1 + ry) / tau;
        b[i] = ry / tau;
        c[i] = ry / tau;
    }
    a[Nx - 1] = 1;
    b[Nx - 1] = 0;
    c[Nx - 1] = 0;

    for (auto i = 0u; i < Nx; i++){
        d[i][0] = 0;
        d[i][Ny - 1] = 0;
    }

    for (auto j = 1u; j < Ny - 1; j++){
        for (auto i = 0u; i < Nx; i++) {
            d[i][j] = 2 * f[i][j] / tau;
        }
    }

    for (auto i = 0u; i < Nx; i++){
        alpha[i][0] = b[0] / a[0];
        betta[i][0] = d[i][0] / a[0];
        for (auto j = 1u; j < Ny; j++) {
            alpha[i][j] = b[j] / (- c[j] * alpha[i][j - 1] + a[j]);
            betta[i][j] = (d[i][j] + c[j] * betta[i][j - 1]) / (- c[j] * alpha[i][j - 1] + a[j]);
        }
    }
    for (auto i = Nx - 1; i >= 0; i--){
        f1[i][Ny - 1] = betta[i][Ny - 1];
        for (auto j = Ny - 2; j >= 0; j--){
            f1[i][j] = alpha[i][j] * f1[i][j + 1] + betta[i][j];
        }
    }
}



void routine(std::vector<std::vector<double>> &f){
    std::vector<std::vector<double>> f1(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> f2(Nx, std::vector<double>(Ny));
    std::ofstream out;
    out.open("f100.txt");
    run_through_x(f, f1, 0);
    run_through_y(f1, f2, 0);

    for (auto k = 0u; k < Nx; k++){
        for (auto l = 0u; l < Ny; l++){
            out << f2[k][l] << ",";
        }
        out << std::endl;
    }

    for (auto i = 1u; i < T / tau; i++){
        run_through_x(f2, f1, i * tau);
        run_through_y(f1, f2, i * tau);
        if (out.is_open()) {
            for (auto k = 0u; k < Nx; k++){
                for (auto l = 0u; l < Ny; l++){
                    out << f2[k][l] << ",";
                }
                out << std::endl;
            }
        }
    }
    out.close();
}

int main(){
    std::vector<std::vector<double>> f(Nx, std::vector<double>(Ny));

    for(auto i = 0u; i < Nx; i++){
        for (auto j = 0u; j < Ny; j++){
            f[i][j] = f_t0(i * hx, j * hy);
        }
    }

    routine(f);

    return 0;
}
