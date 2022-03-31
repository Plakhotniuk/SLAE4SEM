//
// Created by Арсений Плахотнюк on 24.03.2022.
//

#ifndef MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM4_HPP
#define MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM4_HPP

#include <my_project/utility/Overloads.hpp>
#include "../sparse/CSR.hpp"
#include <sstream>
#include <my_project/SlaeBaseException.hpp>
#include <functional>
#include <my_project/matrix/FiveDiagonalMatrix.hpp>
#include "cmath"
#include "my_project/solvers/FiveDiagonalSolver.hpp"
#include "ostream"


double GetСoef1_1(std::function<double(double)>& a, double h, double x){
    return 11/(12*pow(h, 2)) - 0.25*a(x)/h;
}

double GetСoef2_1(std::function<double(double)>& a, std::function<double(double)>& b, double h, double x){
    return -5 / (3*pow(h, 2)) - 5/(6*h) * a(x) + b(x);
}

double GetСoef3_1(std::function<double(double)>& a, double h, double x){
    return 0.5/pow(h, 2) - 1.5 * a(x) / h;
}

double GetCoef4_1(std::function<double(double)>& a, double h, double x){
    return 1/(3*pow(h, 2)) - 0.5/h * a(x);
}

double GetCoef5_1(std::function<double(double)>& a, double h, double x){
    return 1/(12*pow(h, 2)) - 1/(12*h) * a(x);
}

double GetСoef1_i(std::function<double(double)>& a, double h, double x){
    return -1/(12*pow(h, 2)) + a(x)/(12*h);
}

double GetСoef2_i(std::function<double(double)>& a, double h, double x){
    return 4/(3 * pow(h, 2)) - 2/(3*h) * a(x);
}

double GetСoef3_i(std::function<double(double)>& a,std::function<double(double)>& b, double h, double x){
    return b(x) - 2.5/pow(h, 2);
}

double GetCoef4_i(std::function<double(double)>& a, double h, double x){
    return 4/(3*pow(h, 2)) + 2/(3*h) * a(x);
}

double GetCoef5_i(std::function<double(double)>& a, double h, double x){
    return -1/(12*pow(h, 2)) - 1/(12*h) * a(x);
}

double GetСoef1_n_2(std::function<double(double)>& a, double h, double x){
    return 1/(12*pow(h, 2)) - a(x)*1/(12*h);
}

double GetСoef2_n_2(std::function<double(double)>& a, double h, double x){
    return 1 / (3*pow(h, 2)) - 1/(2*h) * a(x);
}

double GetСoef3_n_2(std::function<double(double)>& a, double h, double x){
    return 0.5/pow(h, 2) - 1.5 * a(x) / h;
}

double GetCoef4_n_2(std::function<double(double)>& a, std::function<double(double)>& b, double h, double x){
    return -5/(3*pow(h, 2)) + 5/(6*h) * a(x) + b(x);
}

double GetCoef5_n_2(std::function<double(double)>& a, double h, double x){
    return 11/(12*pow(h, 2)) - 1/(4*h) * a(x);
}

std::vector<double> NonlinearBoundaryValueProblem4(double left_bound_x,double right_bound_x, double left_bound_y,
                                                  double right_bound_y, int number_of_splits,
                                                  std::function<double(double)>& a, std::function<double(double)>& b,
                                                  std::function<double(double)>& f) {
    auto h = (right_bound_x - left_bound_x) / number_of_splits;
    std::vector<double> x(number_of_splits + 1);

    for(int i = 0; i < x.size(); ++i){
        x[i] = left_bound_x + h * i;
    }

    Slae::Matrix::FiveDiagonalMatrix data = Slae::Matrix::FiveDiagonalMatrix(number_of_splits + 1);
    data(0, 1) = 1;
    data(data.rows() - 1, 1) = 1;
    data.fill_row(1, GetСoef1_1(a, h, x[1]), GetСoef2_1(a, b, h, x[1]),
                  GetСoef3_1(a, h, x[1]),GetCoef4_1(a, h, x[1]), GetCoef5_1(a, h, x[1]));

    data.fill_row(data.rows() - 2, GetСoef1_n_2(a, h, x[data.rows() - 2]),
                  GetСoef2_n_2(a, h, x[data.rows() - 2]), GetСoef3_n_2(a, h, x[data.rows() - 2]),
                  GetCoef4_n_2(a, b, h, x[data.rows() - 2]), GetCoef5_n_2(a, h, x[data.rows() - 2]));


    for(int i = 2; i < data.rows() - 2; ++i){
        data.fill_row(i, GetСoef1_i(a, h, x[i]), GetСoef2_i(a, h, x[i]), GetСoef3_i(a, b, h, x[i]),
                      GetCoef4_i(a, h, x[i]), GetCoef5_i(a, h, x[i]));
    }

    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = f(x[i]);
    }
    y[y.size() - 1] = right_bound_y;

    data.multiply_row_by_value(1, data(2, 4)/data(1, 4));
    y[1] *= data(2, 4)/data(1, 4);
    data.substract_row1_from_row2(1, 2);
    y[1] -= y[2];
    data.fill_row(1, 0, data(1, 0), data(1, 1), data(1, 2), data(1, 3));

    data.multiply_row_by_value(data.rows() - 2, data(data.rows() - 3, 0)/data(data.rows() - 2, 0));
    y[data.rows() - 2] *= data(data.rows() - 3, 0)/data(data.rows() - 2, 0);
    data.substract_row1_from_row2(data.rows() - 2, data.rows() - 3);
    y[data.rows() - 2] -= y[data.rows() - 3];

    data.fill_row(data.rows() - 2, data(data.rows() - 2, 1),
                  data(data.rows() - 2, 2), data(data.rows() - 2, 3),
                  data(data.rows() - 2, 4), 0);


    return Slae::Solvers::solveFiveDiagonal(data, y);
}

#endif //MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEM4_HPP
