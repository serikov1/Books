#include <iostream>
#include "cmath"
#include "cstdint"
#include "vector"

double L = 1;
uint8_t N = 100;
double h = 0.01;
double tau = 0.01;
uint16_t T = 2;
uint16_t Nt = 200;
double x0 = 0;

std::vector<double> compute_corner(std::vector<double>& y0, std::vector<double>& y)
{
    for(uint8_t i = 0; i < N-1; i++)
    {
        y[i] = y0[i] - tau * (y0[i + 1] - y0[i]) / h;
    }
    y[N-1] = y[N-2];
    return y;
}


std::vector<std::vector<double>> corner()
{
    std::vector<double> y0(N);
    std::vector<std::vector<double>> y(N, std::vector<double>(N));
    double x = x0;

    for(uint8_t i = 0; i < N/2; i++)
    {
        y0[i] = 1 ;
        x += h;
    }
    for(uint8_t i = N/2; i < N; i++)
    {
        y0[i] = 0.5 ;
        x += h;
    }

    compute_corner(y0, y[0]);
    for(int i = 0; i < N-1; i++)
    {
        compute_corner(y[i], y[i+1]);
    }

    return y;
}


int main() {
    std::vector<std::vector<double>> solution_corner = corner();

    for(auto n = 0; n < N; n++)
    {
        for(uint8_t i = 0; i < N; i++)
        {
            std::cout<< solution_corner[n][i]<< ", " <<std::endl;
        }
    }
    return 0;
}

