//
// Created by Арсений Плахотнюк on 11.03.2022.
//

#ifndef MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM2_HPP
#define MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM2_HPP
#include <my_project/utility/Overloads.hpp>
#include "../sparse/CSR.hpp"
#include <sstream>
#include <my_project/SlaeBaseException.hpp>
#include <functional>
#include <my_project/matrix/ThreeDiagonalMatrix.hpp>
#include "cmath"
#include "my_project/solvers/ThreeDiadonalSolver.hpp"
#include "ostream"

double GetCoef1(std::function<double(double)>& a, double h, double x){
    return 1/(h*h) - a(x)/(2*h);
}


double GetCoef2(std::function<double(double)>& b, double h, double x){
    return -2/(h*h) + b(x);
}


double GetCoef3(std::function<double(double)>& a, double h, double x){
    return 1/(h*h) + a(x)/(2*h);
}

std::vector<double> NonlinearBoundaryValueProblem2(double left_bound_x, double right_bound_x, double left_bound_y,
                                                   double right_bound_y, int number_of_splits,
                                                   std::function<double(double)>& a, std::function<double(double)>& b,
                                                   std::function<double(double)>& f) {
    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;
    std::vector<double> x(number_of_splits + 1);
    for(int i = 0; i < x.size(); ++i){
        x[i] = left_bound_x + h * i;
    }
    // Матрица коэффициентов
    Slae::Matrix::ThreeDiagonalMatrix data = Slae::Matrix::ThreeDiagonalMatrix(number_of_splits + 1);
    // Коэф в начале
    data(0, 1) = 1;
    //Коэф в конце
    data(data.rows() - 1, 1) = 1;
    //Заполение остальных
    for(int i = 1; i < data.rows() - 1; ++i){
            data.fill_row(i, GetCoef1(a, h, x[i]), GetCoef2(b, h, x[i]), GetCoef3(a, h, x[i]));
    }
    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = f(x[i]);
    }
    y[y.size() - 1] = right_bound_y;
    return Slae::Solvers::solveThreeDiagonal(data, y);
}
#endif //MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM2_HPP
