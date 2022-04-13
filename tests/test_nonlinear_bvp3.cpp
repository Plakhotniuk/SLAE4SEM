//
// Created by Арсений Плахотнюк on 14.04.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include "my_project/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>
#include <algorithm>

TEST(NONLINEARBOUNDARYVALUEPROBLEM, TEST) {
std::function<double(double)> f = [](double x) { return 0.; };

std::function<double(double, double, double)> a = [](double x, double y, double y_dev) { return y; };

std::function<double(double, double, double)> b = [](double x, double y, double y_dev) { return 0.; };

std::function<double(double)> y_func = [](double x) { return 2 * tanh(x); };

double left_bound_x = 0.;
double right_bound_x = 1.;
double left_bound_y = 0.;
double right_bound_y = 2. * tanh(1.);
int n = 10;

int max_number_of_splits = 100;
std::vector<double> y(max_number_of_splits + 1);
std::fstream file;
file.open("test_nonlin_3.txt", std::fstream::out);
double h = (right_bound_x - left_bound_x)/max_number_of_splits;
for(int i = 0; i < y.size(); ++i){
y[i] = y_func(left_bound_x + h * i);
}

std::vector<std::vector<double>> solution(n);

std::vector<double> y_0 = InitialApproach(left_bound_x, right_bound_x,
                    left_bound_y, right_bound_y,
                    max_number_of_splits);

//x
for(int j=0; j < y.size(); ++j){
    file << left_bound_x + h * j << " ";
}
file << '\n';

//analitic solution
for(int j=0; j < y.size(); ++j){
    file << y[j] << " ";
}
file << '\n';

//approach
solution[0] = y_0;
for(int j=0; j < y.size(); ++j){
        file << solution[0][j] << " ";
}
file << '\n';
//<< left_bound_x + h * j << " "<<y[j];
for(int i = 1; i < n; ++i){
    std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
            ExpMatrixNonLinearBVP3(left_bound_x, right_bound_x,
                                   left_bound_y, right_bound_y,
                                   max_number_of_splits, a, b, f, solution[i-1]);
    solution[i] = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);
    for(int j=0; j < y.size(); ++j){
        file << solution[i][j] << " ";
    }
    file << '\n';
}

file.close();
}
