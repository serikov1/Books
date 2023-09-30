#include <iostream>
#include <cmath>

// задача 2. Вычислить интеграл I = int|0 - 3| { sin(100x) exp(-x^2) cos(2x) dx}
// использую формулу интегрирования Гаусса Sigma|i=1 - n| { C_i f(x_i) }, где n = 3
const int n = 3;

// узлы и веса для полинома Лежандра
const double Xi[n]={-0.7745967,0,0.7745967};
const double Ci[n]={0.5555556,0.8888889,0.5555556};

// исходная функция
double f(double x) {
    return sin(100 * x) * exp(-pow(x, 2)) * cos(2 * x);
}

// составляем функцию Гаусса со смещением координат из (a, b) в (-1, 1)
double GaussFunc(double a, double b) {
    double diff=(b-a)/2;
    double sum=(a+b)/2;
    double newArgument, newIntegral = 0.0;
    for(uint8_t i = 0; i < n ; i++) {
        newArgument = sum + diff * Xi[i];
        newIntegral += Ci[i] * f(newArgument);
    }
    return diff * newIntegral;
}

int main() {
    double a = 0;
    double b = 3;
    double answer = 0.0;

    // разбиение для интегрирования
    const uint32_t N = 1000000;
    for(uint32_t i = 0; i < N; i++) {
        answer += GaussFunc(a + i * (b-a)/N, a + (i+1) * (b-a)/N);
    }

    printf("Answer is: %.16f", answer);

    return 0;
}
