#include <iostream>
#include <cmath>
#include <vector>

double f1(double xn, double yn) {
    return pow(xn, 2) + pow(yn,2) - 1;
}

double f2(double xn, double yn) {
    return yn - tan(xn);
}

double det(double x, double y) {
    return  1/( 2*x + 2*y/(pow(cos(x), 2)) );
}

std::vector<double> JF(std::vector<double> Xn) {
    return {det(Xn[0], Xn[1])*( f1(Xn[0], Xn[1]) - 2*Xn[1]*f2(Xn[0], Xn[1]) ),
            det(Xn[0], Xn[1])*( f1(Xn[0], Xn[1])/ pow(cos(Xn[0]), 2) + 2*Xn[1]*f2(Xn[0], Xn[1])) };
}


int main() {
    std::vector<double> Xn = {-0.7, -0.7};
    std::vector<double> Xn1 = {0, 0};

    double epsilon = 0.00000001;
    uint32_t n = 1;
    while ( sqrt( std::abs(pow(Xn1[0],2) + pow(Xn1[1], 2)  -  pow(Xn[0],2) - pow(Xn[1], 2) ) ) > epsilon ) {
        Xn1[0] = Xn[0] - JF(Xn)[0];
        Xn1[1] = Xn[1] - JF(Xn)[1];
        if(n > 1) {
            Xn[0] = Xn1[0];
            Xn[1] = Xn1[1];
        }
        n++;
    }

    std::cout<<"Первый нуль функции: " << "x = " << Xn1[0]<< ";" << " y = " << Xn1[1] << "\n" << "Получен на " << n << " шаге МПИ \n" << std::endl;


    Xn = {0.7, 0.7};
    Xn1 = {0, 0};

    n = 1;
    while (sqrt( std::abs(pow(Xn1[0],2) + pow(Xn1[1], 2)  -  pow(Xn[0],2) - pow(Xn[1], 2) ) ) > epsilon ) {
        Xn1[0] = Xn[0] - JF(Xn)[0];
        Xn1[1] = Xn[1] - JF(Xn)[1];
        if(n > 1) {
            Xn[0] = Xn1[0];
            Xn[1] = Xn1[1];
        }
        n++;
    }

    std::cout<<"Второй нуль функции: " << "x = " << Xn1[0]<< ";" << " y = " << Xn1[1] << "\n" << "Получен на " << n << " шаге МПИ \n" << std::endl;

}
