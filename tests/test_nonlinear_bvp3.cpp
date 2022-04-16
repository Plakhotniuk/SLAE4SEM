//
// Created by Арсений Плахотнюк on 14.04.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include "my_project/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include "my_project/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>


TEST(NONLINEARBOUNDARYVALUEPROBLEM1, NONLINEARBOUNDARYVALUEPROBLEM1) {
std::function<double(double, double)> f = [](double x, double y) { return 0.; };

std::function<double(double, double, double)> a = [](double x, double y, double y_dev) { return y; };

std::function<double(double, double, double)> b = [](double x, double y, double y_dev) { return 0.; };

std::function<double(double)> y_func = [](double x) { return 2 * tanh(x); };

double left_bound_x = 0.;
double right_bound_x = 1.;
double left_bound_y = 0.;
double right_bound_y = 2. * tanh(1.);

int n_iterations = 7;
int max_number_of_splits = 100;
std::vector<double> y(max_number_of_splits + 1);
double h = (right_bound_x - left_bound_x)/max_number_of_splits;

std::fstream file;
file.open("test_nonlin_3.txt", std::fstream::out);

//x
//for(int j=0; j < y.size(); ++j){
//    file << left_bound_x + h * j << " ";
//}
//file << '\n';
//
//analytical solution
for(int i = 0; i < y.size(); ++i){
y[i] = y_func(left_bound_x + h * i);
file << y[i] << " ";
}
file << '\n';

//for(int k = 10; k < max_number_of_splits; k+=10){
    std::pair<std::vector<double>, std::vector<double>> solution;

    std::pair<std::vector<double>, std::vector<double>> y_0_y_0_der = InitialApproach_y(left_bound_x, right_bound_x,
                                                                                        left_bound_y, right_bound_y,
                                                                                        max_number_of_splits, y_func);

//approach
    solution = y_0_y_0_der;

    for(int i = 1; i < n_iterations; ++i){
        std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
                ExpMatrixNonLinearBVP3(left_bound_x, right_bound_x,
                                       left_bound_y, right_bound_y,
                                       max_number_of_splits, a, b, f, solution);
        solution.first = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);
    }
    for(int i = 0; i < solution.first.size(); ++i){
        file<< solution.first[i]<< " ";

    }
    file<<"\n";
    for(int i = 0; i < solution.first.size(); ++i){

        file << left_bound_x + h * i << " ";
    }
    file<<"\n";
//}


file.close();
}

TEST(NONLINEARBOUNDARYVALUEPROBLEM_ERROR2, NONLINEARBOUNDARYVALUEPROBLEM_ERROR2) {

    std::function<double(double, double, double)> a = [](double x, double y, double y_der) { return y_der / (2*y); };

    std::function<double(double, double, double)> b = [](double x, double y, double y_der) { return 0.; };

    std::function<double(double, double)> f = [](double x, double y) { return -1/(2*y); };

    std::function<double(double)> y_sol = [](double x) { return log(1.+ x) + 1; };

    double left_bound_x = M_PI / 4 + sqrt(2) / 2;
    double right_bound_x = M_PI;
    double left_bound_y = 1 - sqrt(2) / 2;
    double right_bound_y = 2;

    int n_iterations = 20;
    int max_number_of_splits = 1050;

    std::fstream file;
    file.open("test_nonlin_3_2.txt", std::fstream::out);


for(int k = 20; k < max_number_of_splits; k+=20) {

    //approach
    std::pair<std::vector<double>, std::vector<double>> solution = InitialApproach_y(left_bound_x, right_bound_x,
                                                                                     left_bound_y, right_bound_y,
                                                                                     k, y_sol);;
    double h = (right_bound_x - left_bound_x)/k;
    for (int i = 0; i < n_iterations; ++i)
    {
        std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
                ExpMatrixNonLinearBVP3(left_bound_x, right_bound_x,
                                       left_bound_y, right_bound_y,
                                       k, a, b, f, solution);
        solution.first = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);
        solution.second = Approach_y_der(solution.first, h);
    }


    std::vector<double> err(k + 1);
    for(int i=0; i<k; i++){
        err[i] = abs(solution.first[i] - y_sol(left_bound_x + h * i));
    }
    double max = *std::max_element(err.begin(), err.end());
    file << k << " "<< max << " ";
    file << std::endl;
}
    file.close();
}

TEST(NONLINEARBOUNDARYVALUEPROBLEM3, NONLINEARBOUNDARYVALUEPROBLEM3) {
    std::function<double(double, double, double)> a = [](double x, double y, double y_der) { return -y_der / (2*y); };

    std::function<double(double, double, double)> b = [](double x, double y, double y_der) { return 0.; };

    std::function<double(double, double)> f = [](double x, double y) { return 1/(2*y); };

    std::function<double(double)> y_sol = [](double x) { return x*x; };

    double left_bound_x = 0;
    double right_bound_x = 6;
    double left_bound_y = 5;
    double right_bound_y = 2;

    int n_iterations = 35;
    int max_number_of_splits = 200;
    std::vector<double> y(max_number_of_splits + 1);
    double h = (right_bound_x - left_bound_x)/max_number_of_splits;

    std::fstream file;
    file.open("test_nonlin_3_3.txt", std::fstream::out);

//x
for(int j=0; j < y.size(); ++j){
    file << left_bound_x + h * j << " ";
}
file << '\n';

//for(int k = 10; k < max_number_of_splits; k+=10){

//approach
//{y, y'}
std::pair<std::vector<double>, std::vector<double>> solution = InitialApproach_y(left_bound_x, right_bound_x,
                                                                               left_bound_y, right_bound_y,
                                                                               max_number_of_splits, y_sol);

for(int i = 0; i < n_iterations; ++i){
    std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
            ExpMatrixNonLinearBVP3(left_bound_x, right_bound_x,
                                   left_bound_y, right_bound_y,
                                   max_number_of_splits, a, b, f, solution);
    solution.first = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);

    solution.second = Approach_y_der(solution.first, h);

}

for(int i = 0; i < solution.first.size(); ++i){
    file<< solution.first[i]<< " ";
}
file<<"\n";
    for(int i = 0; i < solution.first.size(); ++i){
        file<< y_sol(left_bound_x + h * i)<< " ";
    }
    file<<"\n";

file.close();
}
