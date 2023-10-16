#include <iostream>
#include <cmath>
#include <vector>

double f1(double x) {
    return pow(x, 2)  - 4 - sin(x);
}

double dev(double x) {
    return  -cos(x) +2*x;
}


int main() {
    double Xn = -2;
    double Xn1 = 0;

    double epsilon = 0.000001;
    uint32_t n = 1;
    while ( Xn1 - Xn > epsilon ) {
        if(n > 1) Xn = Xn1;
        Xn1 = Xn - f1(Xn)/dev(Xn);

        n++;
    }

    std::cout<<"нуль функции: " << "x = " << Xn1 << "\n" << "Получен на " << n << " шаге МПИ \n" << std::endl;

}
