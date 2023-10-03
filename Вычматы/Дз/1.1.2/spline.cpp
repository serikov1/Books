#include "iostream"
#include "vector"
#include "cmath"
#include "iomanip"

// условие задачи
std::vector<double> population = {92228496, 106021537,123202624,132164569,151325798,179323175,203211926,226545805,248709873,281421906};
std::vector<double> years = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};

void Run(std::vector<double> &x, std::vector<double> &f, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, uint32_t N){
    std::vector<double> delta(N + 1);
    std::vector<double> lamda(N + 1);
    std::vector<double> h(N + 1);
    std::vector<double> l(N + 1);

    for (uint32_t k = 1; k <= N; k++) {
        h[k] = x[k] - x[k - 1];
        l[k] = (f[k] - f[k - 1]) / h[k];
    }

    delta[1] = - h[2] / (2 * (h[1] + h[2]));
    lamda[1] = 3 * (l[2] - l[1]) / (2 * (h[1] + h[2]));
    c[0] = 0;
    c[N] = 0;

    for (uint32_t k = 3; k <= N; k++) {
        delta[k - 1] = ( - h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]));
        lamda[k - 1] = ((3 * l[k] - 3 * l[k - 1] - h[k - 1] * lamda[k - 2]) / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]));
    }

    for (uint32_t k = N; k >= 2; k--) {
        c[k - 1] = delta[k - 1] * c [k] + lamda[k - 1];
    }

    for (uint32_t k = 1; k <= N; k++) {
        b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1])  / 3;
        d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
    }
}

double spline (std::vector<double> &x, std::vector<double> &f, double x0) {
    uint32_t N = f.size() - 1;
    std::vector<double> b(N + 1);
    std::vector<double> c(N + 1);
    std::vector<double> d(N + 1);
    Run(x, f, b, c, d, N);

    uint32_t k = 0;
    for (; k < N + 1; k++) {
        if (x0 >= x[k]) continue;
        else break;
    }

    if (x0 >= x[N]) k = N;
    return f[k] + b[k] * (x0 - x[k]) + c[k] * std::pow(x0 - x[k], 2) + d[k] * std::pow(x0 - x[k], 3);
}

int main() {
    uint32_t amount = 308745538;

    std::cout<< std::fixed << std::setprecision(0) << "Экстраполяция сплайнами. Населения Омэрики на 2010г: " << spline(years, population, 2010) <<std::endl;
    std::cout<< std::fixed << std::setprecision(8) << "Отношение экстраполированной величины к реальной: " << spline(years, population, 2010) / amount <<std::endl;
}