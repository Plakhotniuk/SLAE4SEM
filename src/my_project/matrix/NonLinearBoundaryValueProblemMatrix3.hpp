//
// Created by Арсений Плахотнюк on 13.04.2022.
//

#ifndef MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
#define MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
#include "/Users/arseniy/Desktop/SLAE4SEM/src/my_project/utility/Overloads.hpp"
#include "/Users/arseniy/Desktop/SLAE4SEM/src/my_project/sparse/CSR.hpp"
#include <sstream>
#include "/Users/arseniy/Desktop/SLAE4SEM/src/my_project/SlaeBaseException.hpp"
#include <functional>
#include "ThreeDiagonalMatrix.hpp"
#include "cmath"
#include "/Users/arseniy/Desktop/SLAE4SEM/src/my_project/solvers/ThreeDiadonalSolver.hpp"
#include "ostream"

enum COLUMNINDEX {
    C_FIRST = 1,
    C_SECOND = 2,
    C_THIRD = 3,
};

template<int C_Index>
struct Calc{
    static double calc(std::function<double(double, double, double)>& a,
                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der);
};

template<>
double Calc<COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                        std::function<double(double, double, double)>& b,
                                        double h, double x, double y, double y_der){
    return 1/(h*h) - a(x, y, y_der)/(2*h);
}

template<>
double Calc<COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                         std::function<double(double, double, double)>& b,
                                         double h, double x, double y, double y_der){
    return -2/(h*h) + b(x, y, y_der);
}

template<>
double Calc<COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                        std::function<double(double, double, double)>& b,
                                        double h, double x, double y, double y_der){
    return 1/(h*h) + a(x, y, y_der)/(2*h);
}

std::pair<std::vector<double>, std::vector<double>> InitialApproach(double left_bound_x, double right_bound_x,
                                    double left_bound_y, double right_bound_y,
                                    unsigned int number_of_splits){
    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // начальное приближение
    std::vector<double> y_0(number_of_splits + 1);
    for(unsigned int i = 0; i < y_0.size(); ++i)
    {
        y_0[i] = left_bound_y + (left_bound_x + h * i - left_bound_x) /
                                (right_bound_x - left_bound_x) * (right_bound_y - left_bound_y);
    }
    std::vector<double> y_0_der(number_of_splits + 1, (right_bound_y - left_bound_y) / (right_bound_x - left_bound_x));
    return {y_0, y_0_der};
}

std::vector<double> y_der(std::vector<double> y_prev, double h){
    std::vector<double> y_der(y_prev.size());
    for(int i = 0; i < y_prev.size() - 1; ++i){
        y_der[i] = (y_prev[i+1] - y_prev[i]) / h;
    }
    y_der.back() = y_der[y_der.size()-2];
    return y_der;
}

std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>>
ExpMatrixNonLinearBVP3(double left_bound_x, double right_bound_x,
                       double left_bound_y, double right_bound_y,
                       unsigned int number_of_splits,
                       std::function<double(double, double, double)>& a,
                       std::function<double(double, double, double)>& b,
                       std::function<double(double, double)>& f, std::pair<std::vector<double>, std::vector<double>> y_prev){

    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // Матрица коэффициентов
    Slae::Matrix::ThreeDiagonalMatrix data = Slae::Matrix::ThreeDiagonalMatrix(number_of_splits + 1);
    // Коэф в начале
    data(0, 1) = 1;
    //Коэф в конце
    data(data.rows() - 1, 1) = 1;

    //Заполение остальных
    for(unsigned int i = 1; i < data.rows() - 1; ++i){
        data.fill_row(i, Calc<COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]));
    }

    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(unsigned int i = 1; i < y.size() - 1; ++i){
        y[i] = f(left_bound_x + h * i, y_prev.first[i]);
    }
    y.back() = right_bound_y;
    return {data, y};
}

#endif //MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
