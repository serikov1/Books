#include <iostream>
#include "cmath"
#include "cstdint"
#include "vector"

double L = 20.5;
uint8_t N = 42;
double h = L / (N - 1);
double tau = 10*h;
double T = 18;
double x0 = 0;
double a = 1;

double u0(double x)
{
    return sin(4* M_PI * x / L);
}

std::vector<double> compute_corner(std::vector<double>& y0, std::vector<double>& y)
{
    for(uint8_t i = 1; i < N; i++)
    {
        y[i] = y0[i] + tau * (y0[i - 1] - y0[i]) / h;
    }
    y[0] = y[N-1];
    return y;
}

std::vector<double> compute_Lax(std::vector<double>& y0, std::vector<double>& y){
    for (auto i = 1; i < N - 1; i++){
        y[i] = y0[i] + tau * (a * (y0[i - 1] - y0[i + 1]) / (2 * h) + pow(a, 2) * 0.5 * tau * (y0[i + 1] - 2 * y0[i] + y0[i - 1]) / pow(h, 2));
    }
    y[N - 1] = y0[N - 1] + tau * (a * (y0[N - 2] - y0[0]) / (2 * h) + pow(a, 2) * 0.5 * tau * (y0[0] - 2 * y0[N - 1] + y0[N - 2]) / pow(h, 2));
    y[0] = y0[0] + tau * (a * (y0[N - 1] - y0[1]) / (2 * h) + pow(a, 2) * 0.5 * tau * (y0[1] - 2 * y0[0] + y0[N - 1]) / pow(h, 2));
    return y;
}

std::vector<std::vector<double>> corner()
{
    std::vector<double> y0(N);
    std::vector<std::vector<double>> y(N, std::vector<double>(N));
    double x = x0;

    for(uint8_t i = 0; i < N; i++)
    {
        y0[i] = u0(x);
        x += h;
    }

    compute_corner(y0, y[0]);
    for(uint8_t i = 0; i < T-1; i++)
    {
        compute_corner(y[i], y[i+1]);
    }

    return y;
}

std::vector<std::vector<double>> Lax()
{
    std::vector<double> y0(N);
    std::vector<std::vector<double>> y(N, std::vector<double>(N));
    double x = x0;

    for(uint8_t i = 0; i < N; i++)
    {
        y0[i] = u0(x);
        x += h;
    }

    compute_Lax(y0, y[0]);
    for(uint8_t i = 0; i < T-1; i++)
    {
        compute_Lax(y[i], y[i+1]);
    }

    return y;
}

int main() {
    std::vector<std::vector<double>> solution_corner = corner();
    std::vector<std::vector<double>> solution_Lax = Lax();

    for(auto n = 0; n < 1; n++)
    {
//        std::cout<<"t = "<< n << std::endl;
        for(uint8_t i = 0; i < N; i++)
        {
            std::cout<< solution_corner[n][i]<< ", " <<std::endl;
        }
    }
//    for(auto n = 16; n < 17; n++)
//    {
////        std::cout<<"t = "<< n << std::endl;
//        for(uint8_t i = 0; i < N; i++)
//        {
//            std::cout<< solution_Lax[n][i]<< ", " <<std::endl;
//        }
//    }
//    return 0;
}
