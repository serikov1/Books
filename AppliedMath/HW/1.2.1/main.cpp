#include <iostream>
#include <cmath>
#include <vector>
#include <cstdint>

//Построить алгоритм метода пристрелки для вычисления
//решения следующей нелинейной задачи:

// y'' - x * sqrt (y) = 0, 0<x<1
// y(0) = 0, y(1) = 2

double y_1 = 2;
double epsilon = 0.001;

// y' = z = g(x, y, z)
double g(double x, double y, double z) {
    return z;
}

// y'' = f(x, y, z) = g'
double f(double x, double y, double z) {
    return x * std::sqrt(y);
}


class RungeKutte {
public:
    double K1(std::vector<double>& u) const {
        return f(u[0], u[1], u[2]);
    }
    double K2(std::vector<double>& u, double k1, double p1) const {
        return f(u[0] + 0.5 * h, u[1] + 0.5 * h * p1, u[2] + 0.5 * h * k1);
    }
    double K3(std::vector<double>& u, double k2, double p2) const {
        return f(u[0] + 0.5 * h, u[1] + 0.5 * h * p2, u[2] + 0.5 * h * k2);
    }
    double K4(std::vector<double>& u, double  k3, double p3) const {
        return f(u[0] + h, u[1] + h * p3, u[2] + h * k3);
    }

    double P1(std::vector<double>& u) const {
        return g(u[0], u[1], u[2]);
    }
    double P2(std::vector<double>& u, double k1, double p1) const {
        return g(u[0] + 0.5 * h, u[1] + 0.5 * h * p1, u[2] + 0.5 * h * k1);
    }
    double P3(std::vector<double>& u, double k2, double p2) const {
        return  g(u[0] + 0.5 * h, u[1] + 0.5 * h * p2, u[2] + 0.5 * h * k2);
    }
    double P4(std::vector<double>& u, double  k3, double p3) const {
        return g(u[0] + h, u[1] + h * p3, u[2] + h * k3);
    }

    std::vector<double> delta(std::vector<double>& u) {
        double k1 = K1(u);
        double p1 = P1(u);
        double k2 = K2(u, k1, p1);
        double p2 = P2(u, k1, p1);
        double k3 = K3(u, k2, p2);
        double p3 = P3(u, k2, p2);
        double k4 = K4(u, k3, p3);
        double p4 = P4(u, k3, p3);

        return {h, 1./6 * h * (p1 + 2 * p2 + 2 * p3 + p4), 1./6 * h * (k1 + 2 * k2 + 2 * k3 + k4)};
    }

    std::vector<double> calculate(std::vector<double>& u) {
        std::vector<double> v = delta(u);
        u[0] += v[0];
        u[1] += v[1];
        u[2] += v[2];
        return u;
    }

    std::vector<double> getSolution(std::vector<double>& u) {
        std::vector<double> changeable(3);
        changeable = u;
        uint32_t n = 0;
        double z = u[2];
        while (std::abs(changeable[1] - y_1) > epsilon) {
            changeable = {0, 0, z};
            while (h * n < T) {
                changeable = calculate(changeable);
                n++;
            }
            if (changeable[1] > y_1) z -= std::abs(changeable[1] - y_1) ;
            else z += std::abs(changeable[1] - y_1);
            changeable[2] = z;
            n = 0;
        }
//        std::cout<<"z = " << z << std::endl;
        return changeable;
    }

    double geth() const {
        return h;
    }
    double getT() const {
        return T;
    }

    void seth(double h_new) {
        h = h_new;
    }
    void setT(double T_new) {
        T = T_new;
    }

private:
    double h = 0.001;
    double T = 1;
};




int main() {
    std::vector<double> beginCondition = {0, 0, 1};
    RungeKutte ex;
    std::vector<double> res = ex.getSolution(beginCondition);
    std::cout<< "x = " << res[0] << " y = " << res[1]<< " z = " << res[2] <<std::endl;
    return 0;
}
