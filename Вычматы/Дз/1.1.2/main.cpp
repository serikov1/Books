#include <iostream>
#include <cmath>

//Методом простой итерации найти ширину
//функции на полувысоте с точностью 10e-3: x * exp(-x^2), x>=0

double func(double num, double half) {
    return num * exp(-pow(num, 2)) - half;
}

double derivative(double num) {
    return exp(-pow(num, 2)) * (1 - 2 * pow(num, 2));
}



int main() {

    double x0 = 0.6;
    double xn = x0;
    double xn1 = 0;
    double tau = 0.55;
    double epsilon = 0.00001;
    uint32_t n = 1;
    while(std::abs(xn1 - xn) > epsilon) {
        if(n > 1) xn = xn1;
        xn1 = xn + tau * derivative(xn);
        n++;
    }
    
    std::cout<<"Экстремум функции (нуль производной): " << xn1 << "\n" << "Получен на " << n << " шаге МПИ\n" << std::endl;

    double half = func(xn1, 0)/2;
    std::cout<<"Полувысота функции в точке экстремума: " << half << "\n" << std::endl;

    x0 = 0.1;
    tau = -0.5;
    xn = x0;
    xn1 = 0;
    n = 1;
    while(std::abs(xn1 - xn) > epsilon) {
        if(n > 1) xn = xn1;
        xn1 = xn + tau * func(xn, half);
        n++;
    }
    double left = xn1;
    std::cout<<"Левый нуль функции на полувысоте: " << xn1 << "\n" << "Получен на " << n << " шаге МПИ \n" << std::endl;

    x0 = 1.4;
    tau = 0.9;
    xn = x0;
    xn1 = 0;
    n = 1;
    while(std::abs(xn1 - xn) > epsilon) {
        if(n > 1) xn = xn1;
        xn1 = xn + tau * func(xn, half);
        n++;
    }
    std::cout<<"Правый нуль функции на полувысоте: " << xn1 << "\n" << "Получен на " << n << " шаге МПИ \n" << std::endl;
    std::cout<<"Полуширина функци на полувысоте: " << xn1 -  left << "\n" << std::endl;
}
