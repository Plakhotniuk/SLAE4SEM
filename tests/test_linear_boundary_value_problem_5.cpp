//
// Created by Арсений Плахотнюк on 24.03.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include "my_project/matrix/LinearBoundaryValueProblemMatrix5.hpp"
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>
#include <algorithm>


TEST(NONLINEARBOUNDARYVALUEPROBLEM4, TEST) {
    std::function<double(double)> f = [](double x) { return 0.; };

    std::function<double(double)> a = [](double x) { return 0.; };

    std::function<double(double)> b = [](double x) { return 1.; };

    std::function<double(double)> y_func = [](double x) { return cos(x); };

    double left_bound_x = 0.;
    double right_bound_x = M_PI;
    double left_bound_y = 1.;
    double right_bound_y = -1.;

    int max_number_of_splits = 500;
    std::fstream file;
    file.open("test_4_func3.txt", std::fstream::out);
    std::pair<Slae::Matrix::FiveDiagonalMatrix, std::vector<double>> matrix =
            ExpandedMatrixForLinearBoundaryValueProblem5(left_bound_x, right_bound_x,
                                                         left_bound_y, right_bound_y,
                                                         max_number_of_splits, a, b, f);
    std::vector<double> solution = Slae::Solvers::solveFiveDiagonal(matrix.first, matrix.second);

    std::vector<double> y(max_number_of_splits + 1);

    double h = (right_bound_x - left_bound_x)/max_number_of_splits;
    for(int i = 0; i < y.size(); ++i){
        y[i] = y_func(left_bound_x + h * i);
    }
    std::vector<double> err(max_number_of_splits);
    std::cout << "[";
    for(int i=0; i<max_number_of_splits; i++){
        err[i] = abs(solution[i] - y[i]);
        std::cout << solution[i] << " ";
        file << solution[i] << " "<< left_bound_x + h * i <<" "<< y[i];
        file << std::endl;
    }
    std::cout << "]";
    std::cout<< std::endl;
    file << std::endl;
    file.close();
}

TEST(MAXERROR4, TEST){
    std::function<double(double)> f = [](double x) { return 2*x; };

    std::function<double(double)> a = [](double x) { return 0.; };

    std::function<double(double)> b = [](double x) { return -1.; };

    std::function<double(double)> y_func = [](double x) { return sinh(x) / sinh(1.) - 2*x; };

    double left_bound_x = 0.;
    double right_bound_x = 1.;
    double left_bound_y = 0.;
    double right_bound_y = -1.;

    int max_number_of_splits = 220;
    std::fstream file;
    file.open("test_4_err3.txt", std::fstream::out);

    for(int j = 20; j < max_number_of_splits; j += 10){
        std::pair<Slae::Matrix::FiveDiagonalMatrix, std::vector<double>> matrix =
                ExpandedMatrixForLinearBoundaryValueProblem5(left_bound_x, right_bound_x, left_bound_y,
                                                             right_bound_y, j, a, b, f);
        std::vector<double> solution = Slae::Solvers::solveFiveDiagonal(matrix.first, matrix.second);
        std::vector<double> y(j + 1);

        double h = (right_bound_x - left_bound_x)/j;
        for(int i = 0; i < y.size(); ++i){
            y[i] = y_func(left_bound_x + h * i);
        }
        std::vector<double> err(j);
        for(int i=0; i<j; i++){
            err[i] = abs(solution[i] - y[i]);
        }
        double max = *std::max_element(err.begin(), err.end());
        file << j << " "<< max << " ";
        file << std::endl;
    }
    file.close();
}


