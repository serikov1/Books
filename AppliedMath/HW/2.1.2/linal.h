#ifndef INC_2_1_2_LINAL_H
#define INC_2_1_2_LINAL_H

#include "vector"

std::vector<std::vector<double>>& matrix_to_matrix(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& AB);
std::vector<double>& matrix_to_vector(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& Ab);
std::vector<double>& vector_minus_vector(std::vector<double>& a, std::vector<double>& b, std::vector<double>& a_b);
std::vector<double>& vector_plus_vector(std::vector<double>& a, std::vector<double>& b, std::vector<double>& a_b);

#endif //INC_2_1_2_LINAL_H
