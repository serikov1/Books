#include "iostream"
#include "vector"
#include "iomanip"

// По приведенным данным построить интерполянт в форме Ньютона. Вычислить экстраполированное значение численности населения
// США в 2010 году и сравнить с точным значением 308 745 538 человек.
// По этим же данным построить сплайн-аппроксимацию, экстраполировать данные на 2010 год, сравнить с точным значением.

// условие задачи
std::vector<double> population = {92228496, 106021537,123202624,132164569,151325798,179323175,203211926,226545805,248709873,281421906};
std::vector<double> years = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};

// вычисление разделенной разности
double Difference (std::vector<double> &x, std::vector<double> &f, uint16_t k){
    double y_k = 0;
    for (uint32_t i = 0; i <= k; i++) {
        double product = 1;
        for (uint32_t j = 0; j <= k; j++) {
            if (j == i) continue;
            product *= (x[i] - x[j]);
        }
        y_k += (f[i] / product);
    }
    return y_k;
}

double NewtonFunc (std::vector<double> &x, std::vector<double> &f, double x0){
    double P = f[0];
    for (uint32_t k = 1; k < f.size(); k++) {
        double product = 1;
        for (uint32_t i = 0; i < k ; i++) {
            product *= (x0 - x[i]);
        }
        P += (Difference(x, f, k) * product);
    }
    return P;
}

int main() {
    uint32_t amount = 308745538;

    std::cout<< std::fixed << std::setprecision(0) << "Экстраполяция сплайнами. Населения Омэрики на 2010г: " << NewtonFunc(years, population, 2010) <<std::endl ;
    std::cout<< std::fixed << std::setprecision(8) << "Отношение экстраполированной величины к реальной: "<< NewtonFunc(years, population, 2010) / amount << "\n" <<std::endl;
}