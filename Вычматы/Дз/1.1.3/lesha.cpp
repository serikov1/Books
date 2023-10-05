#include "iostream"
#include "vector"
#include "cmath"
#include <fmt/core.h>

std::vector<double> population = {92228496, 106021537,123202624,132164569,151325798,179323175,203211926,226545805,248709873,281421906};
std::vector<double> grid = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};

std::vector<double> y2 = {0, 1, 8, 27};
std::vector<double> x2 = {0, 1, 2, 3};

double difference (std::vector<double> &x, std::vector<double> &f, uint16_t k){
    double y_k = 0;
    for (auto i = 0u; i <= k; i++){
        double product = 1;
        for (auto j = 0u; j <= k; j++){
            if (j == i) continue;
            product *= (x[i] - x[j]);
        }
        y_k += (f[i] / product);
    }
    return y_k;
}


std::vector<double> coefs;
double Newton (std::vector<double> &x, std::vector<double> &f, double x0){
    coefs.resize(f.size() + 1, 0);
    std::vector<double> mid(f.size() + 1, 0);
    mid[0] = -x[0];
    mid[1] = 1;
    uint16_t check = 0;
    double P = f[0];
    double diff;
    for (auto k = 1u; k < f.size(); k++){
        double product = 1;
        diff = difference(x, f, k);
        for (auto i = 0u; i <= k - 1; i++){
            product *= (x0 - x[i]);
            if (check < i) {
                for (auto j = check + 2; j > 0; j--) {
                    mid[j] = mid[j - 1];
                }
                for (auto j = 1u; j < check + 2; j++) {
                    mid[j] -= (x[i] * mid[j + 1]);
                }
                mid[0] *= (-x[i]);

                check++;
            }
        }
        for (auto l = 0; l <= check + 2; l++) {
            coefs[l] += (diff * mid[l]) ;
        }
        P += (diff * product);
    }
    coefs[0] += f[0];
    return P;
}



void run_through(std::vector<double> &x, std::vector<double> &f, std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, uint32_t N){
    std::vector<double> delta(N + 1);
    std::vector<double> lamda(N + 1);
    std::vector<double> h(N + 1);
    std::vector<double> l(N + 1);

    for (auto k = 1u; k <= N; k++){
        h[k] = x[k] - x[k - 1];
        l[k] = (f[k] - f[k - 1]) / h[k];
    }
    delta[1] = - h[2] / (2 * (h[1] + h[2]));
    lamda[1] = 3 * (l[2] - l[1]) / (2 * (h[1] + h[2]));
    c[0] = 0;
    c[N] = 0;
    for (auto k = 3u; k <= N; k++){
        delta[k - 1] = ( - h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]));
        lamda[k - 1] = ((3 * l[k] - 3 * l[k - 1] - h[k - 1] * lamda[k - 2]) / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]));
    }
    for (auto k = N; k >= 2; k--) {
        c[k - 1] = delta[k - 1] * c [k] + lamda[k - 1];
    }
    for (auto k = 1u; k <= N; k++) {
        b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1])  / 3;
        d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
    }
}

double spline (std::vector<double> &x, std::vector<double> &f, double x0){
    uint32_t N = f.size() - 1;
    std::vector<double> a(N + 1);
    std::vector<double> b(N + 1);
    std::vector<double> c(N + 1);
    std::vector<double> d(N + 1);
    run_through(x, f, a, b, c, d, N);

    uint32_t k = 0;
    for (; k < N + 1; k++) {
        if (x0 >= x[k]) continue;
        else break;
    }
    if (x0 >= x[N]) k = N;
    return f[k] + b[k] * (x0 - x[k]) + c[k] * std::pow(x0 - x[k], 2) + d[k] * std::pow(x0 - x[k], 3);
}



int main(){
    double amount = 308745538;

    fmt::print("Метод Ньютона\nP({}) = {} - экстраполированное значение населения на 2010г\n",2010 , Newton(grid, population, 2010));
    fmt::print("P({}) / {} = {} - отношение экстраполированной величины к реальной \n\n",2010, amount, Newton(grid, population, 2010) / amount);

    fmt::print("Сплайны\nF({}) = {} - экстраполированное значение населения на 2010г \n", 2010, spline(grid, population, 2010));
    fmt::print("F({}) / {} = {} - отношение экстраполированной величины к реальной\n", 2010, amount, spline(grid, population, 2010) / amount );

//    fmt::print("Ньютон {}\n", Newton(x2, y2, 9));
//    for (auto j = 0; j < coefs.size(); j++) {
//        fmt::print("c(x^{}) = {}\n", j, coefs[j]);
//    }
}