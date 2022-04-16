//
// Created by Арсений Плахотнюк on 17.04.2022.
//

#ifndef MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX5_HPP
#define MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX5_HPP


#include "my_project/utility/Overloads.hpp"
#include "my_project/sparse/CSR.hpp"
#include <sstream>
#include "my_project/SlaeBaseException.hpp"
#include <functional>
#include "FiveDiagonalMatrix.hpp"
#include "cmath"
#include "my_project/solvers/FiveDiagonalSolver.hpp"
#include "ostream"

enum ROWINDEX {
    R_FIRST = 1,
    R_I = 2,
    R_N_2 = 3,
    R_THREE_DIAG = 4
};

enum COLUMNINDEX {
    C_FIRST = 1,
    C_SECOND = 2,
    C_THIRD = 3,
    C_FOURTH = 4,
    C_FIFTH = 5
};

template<int R_Index, int C_Index>
struct Calc{
    static double calc(std::function<double(double, double, double)>& a,
                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der);
};

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                                           std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 11/(12*(h*h)) - 0.25*a(x, y, y_der)/h;
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                                            std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -5 / (3*(h*h)) - 5/(6*h) * a(x, y, y_der) + b(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                                           std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 0.5/(h*h) + 1.5 * a(x, y, y_der) / h;
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double, double, double)>& a,
                                                            std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 1/(3*(h*h)) - 0.5/h * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double, double, double)>& a,
                                                           std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -1/(12*(h*h)) + 1/(12*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -1/(12*(h*h)) + a(x, y, y_der)/(12*h);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                                        std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 4/(3 * (h*h)) - 2/(3*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return b(x, y, y_der) - 2.5/(h*h);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double, double, double)>& a,
                                                        std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 4/(3*(h*h)) + 2/(3*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double, double, double)>& a,
                                                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -1/(12*(h*h)) - 1/(12*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                                         std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -1/(12*(h*h)) - a(x, y, y_der)*1/(12*h);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                                          std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 1 / (3*(h*h)) + 1/(2*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                                         std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 0.5/(h*h) - 1.5 * a(x, y, y_der) / h;
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double, double, double)>& a,
                                                          std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -5/(3*(h*h)) + 5/(6*h) * a(x, y, y_der) + b(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double, double, double)>& a,
                                                         std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 11/(12*(h*h)) + 1/(4*h) * a(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                                                std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 1/(h*h) - a(x, y, y_der)/(2*h);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                                                 std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return -2/(h*h) + b(x, y, y_der);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                                                std::function<double(double, double, double)>& b, double h, double x, double y, double y_der){
    return 1/(h*h) + a(x, y, y_der)/(2*h);
}

std::vector<double> Approach_y_der5(std::vector<double> y_prev, double h){
    std::vector<double> y_der(y_prev.size());
    y_der[0] = (y_prev[1] - y_prev[0]) / h;
    y_der[1] = (y_prev[2] - y_prev[0]) / (2*h);

    for(int i = 2; i < y_prev.size() - 3; ++i)
        y_der[i] = (-y_prev[i+2] + 8*y_prev[i+1] - 8*y_prev[i-1] + y_prev[i-2]) / (12*h);

    y_der[y_prev.size() - 2] = (y_prev.back() - y_prev[y_prev.size() - 3]) / (2*h);
    y_der.back() = (y_prev.back() - y_prev[y_prev.size() - 2]) / h;

    return y_der;
}

std::pair<std::vector<double>, std::vector<double>> InitialApproach_y5(double left_bound_x, double right_bound_x,
                                                                      double left_bound_y, double right_bound_y,
                                                                      unsigned int number_of_splits,
                                                                      std::function<double(double)>& func){
    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // начальное приближение
    std::vector<double> y_0(number_of_splits + 1);
    for(unsigned int i = 0; i < y_0.size(); ++i)
    {
//        y_0[i] = func(left_bound_x + h * i);
        y_0[i] = left_bound_y + (left_bound_x + h * i - left_bound_x) /
                                (right_bound_x - left_bound_x) * (right_bound_y - left_bound_y);
    }
    std::vector<double> y_0_der(number_of_splits + 1, (right_bound_y - left_bound_y) / (right_bound_x - left_bound_x));
//    std::vector<double> y_0_der = Approach_y_der(y_0, h);
    return {y_0, y_0_der};
}


std::pair<Slae::Matrix::FiveDiagonalMatrix, std::vector<double>>
ExpMatrixNonLinearBVP5(double left_bound_x, double right_bound_x, double left_bound_y,
                       double right_bound_y, int number_of_splits,
                       std::function<double(double, double, double)>& a, std::function<double(double, double, double)>& b,
                       std::function<double(double, double)>& f, std::pair<std::vector<double>, std::vector<double>> y_prev) {
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    Slae::Matrix::FiveDiagonalMatrix data = Slae::Matrix::FiveDiagonalMatrix(number_of_splits + 1);
    data(0, 2) = 1;
    data(data.rows() - 1, 2) = 1;

    //Использование коэффициентов для 3х диагональной матрицы

    data.fill_row(1, 0,
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * 1, y_prev.first[1], y_prev.second[1]),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * 1, y_prev.first[1], y_prev.second[1]),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * 1, y_prev.first[1], y_prev.second[1]), 0);

    data.fill_row(data.rows() - 2, 0,
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * (data.rows() - 2),
                                                                           y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * (data.rows() - 2),
                                                                            y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * (data.rows() - 2),
                                                                           y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]), 0);

    //Элементарные преобразования и заполнение коэффициентов 5ти диагональной матрицы

    // 1 < i < n - 2
    for(int i = 2; i < data.rows() - 2; ++i){
        data.fill_row(i, Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]));
    }
    //i = 1
//    data.fill_row(1, 0, Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h, y_prev.first[1], y_prev.second[1]),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h, y_prev.first[1], y_prev.second[1]),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h, y_prev.first[1], y_prev.second[1]),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h, y_prev.first[1], y_prev.second[1]));
//
//    data.multiply_row_by_value(1, Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 2, y_prev.first[2], y_prev.second[2]) /
//            Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h, y_prev.first[1], y_prev.second[1]));
//    data.fill_row(1, 0, data(1, 1) - data(2, 0),
//                  data(1, 2) - data(2, 1),
//                  data(1, 3) - data(2, 2),
//                  data(1, 4) - data(2, 3));
//
//    // i = n - 2
//    data.fill_row(data.rows() - 2, Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h *(data.rows() - 2),
//    y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]),
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h *(data.rows() - 2),
//                  y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]) ,
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h *(data.rows() - 2),
//                  y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]),
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h *(data.rows() - 2),
//                  y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]), 0);
//
//
//    data.multiply_row_by_value(data.rows() - 2,
//                               Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h *(data.rows() - 3),
//                               y_prev.first[data.rows() - 3], y_prev.second[data.rows() - 3]) /
//                                       Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h *(data.rows() - 2),
//                                       y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]));
//
//    data.fill_row(data.rows() - 2, data(data.rows() - 2, 0) - data(data.rows() - 3, 1),
//                  data(data.rows() - 2, 1) - data(data.rows() - 3, 2),
//                  data(data.rows() - 2, 2) - data(data.rows() - 3, 3),
//                  data(data.rows() - 2, 3) - data(data.rows() - 3, 4), 0);

    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = f(left_bound_x + h * i, y_prev.first[i]);
    }

    y.back() = right_bound_y;

    y[1] *= Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 2, y_prev.first[2], y_prev.second[2]) /
            Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 1, y_prev.first[1], y_prev.second[1]);

    y[1] -= y[2];

    y[data.rows() - 2] *= Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc( a, b, h,left_bound_x + h * (data.rows() - 3), y_prev.first[data.rows() - 3], y_prev.second[data.rows() - 3]) /
                          Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * (data.rows() - 2), y_prev.first[data.rows() - 2], y_prev.second[data.rows() - 2]);

    y[data.rows() - 2] -= y[data.rows() - 3];

    return {data, y};
}


#endif //MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX5_HPP
