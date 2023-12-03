#include <iostream>
#include <cmath>
#include <vector>

double tau = 0.1;
double u_0 = 1;
double t_0 = 0;
double t_1 = 10;
double l = t_1/tau;
std::vector<double> u_i;
std::vector<double> u_i1;
int main() {
    u_i.push_back(u_0);
    for(auto i = 0; i < l; i++) {
        u_i1.push_back(u_i[i]*(tau + 1));
        u_i.push_back(u_i1[i]);
        std::cout<<u_i1[i]<<", ";
    }
    return 0;
}
