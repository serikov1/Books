#include "linal.h"

std::vector<std::vector<double>>& matrix_to_matrix(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& AB)
{
    for(auto i = 0u; i < 3; i++){
        for(auto k = 0u; k < 3; k++) {
            for (auto j = 0u; j < 3; j++) {
                AB[i][k] += A[i][j] * B[j][k];
            }
        }
    }

    return AB;
}

std::vector<double>& matrix_to_vector(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& Ab)
{
    for(auto i = 0u; i < 3; i++){
        for(auto j = 0u; j < 3; j++){
            Ab[i] += A[i][j] * b[j];
        }
    }

    return Ab;
}

std::vector<double>& vector_minus_vector(std::vector<double>& a, std::vector<double>& b, std::vector<double>& a_b){
    for(auto i = 0u; i < 3; i++){
        a_b[i] = a[i] - b[i];
    }

    return a_b;
}

std::vector<double>& vector_plus_vector(std::vector<double>& a, std::vector<double>& b, std::vector<double>& a_b){
    for(auto i = 0u; i < 3; i++){
        a_b[i] = a[i] + b[i];
    }

    return a_b;
}