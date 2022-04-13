//
// Created by Арсений Плахотнюк on 24.03.2022.
//

#ifndef MY_PROJECT_LINEARBOUNDARYVALUEPROBLEM4_HPP
#define MY_PROJECT_LINEARBOUNDARYVALUEPROBLEM4_HPP

#include <my_project/utility/Overloads.hpp>
#include "../sparse/CSR.hpp"
#include <sstream>
#include <my_project/SlaeBaseException.hpp>
#include <functional>
#include <my_project/matrix/FiveDiagonalMatrix.hpp>
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
    static double calc(std::function<double(double)>& a,
                     std::function<double(double)>& b, double h, double x);
};

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIRST>::calc(std::function<double(double)>& a,
                                  std::function<double(double)>& b, double h, double x){
    return 11/(12*(h*h)) - 0.25*a(x)/h;
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_SECOND>::calc(std::function<double(double)>& a,
                                                           std::function<double(double)>& b, double h, double x){
    return -5 / (3*(h*h)) - 5/(6*h) * a(x) + b(x);
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_THIRD>::calc(std::function<double(double)>& a,
                                                            std::function<double(double)>& b, double h, double x){
    return 0.5/(h*h) + 1.5 * a(x) / h;
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double)>& a,
                                                           std::function<double(double)>& b, double h, double x){
    return 1/(3*(h*h)) - 0.5/h * a(x);
}

template<>
double Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double)>& a,
                                                            std::function<double(double)>& b, double h, double x){
    return -1/(12*(h*h)) + 1/(12*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(std::function<double(double)>& a,
                                                           std::function<double(double)>& b, double h, double x){
    return -1/(12*(h*h)) + a(x)/(12*h);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_SECOND>::calc(std::function<double(double)>& a,
                                                       std::function<double(double)>& b, double h, double x){
    return 4/(3 * (h*h)) - 2/(3*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_THIRD>::calc(std::function<double(double)>& a,
                                                        std::function<double(double)>& b, double h, double x){
    return b(x) - 2.5/(h*h);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double)>& a,
                                                       std::function<double(double)>& b, double h, double x){
    return 4/(3*(h*h)) + 2/(3*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double)>& a,
                                                        std::function<double(double)>& b, double h, double x){
    return -1/(12*(h*h)) - 1/(12*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIRST>::calc(std::function<double(double)>& a,
                                                       std::function<double(double)>& b, double h, double x){
    return -1/(12*(h*h)) - a(x)*1/(12*h);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_SECOND>::calc(std::function<double(double)>& a,
                                                        std::function<double(double)>& b, double h, double x){
    return 1 / (3*(h*h)) + 1/(2*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_THIRD>::calc(std::function<double(double)>& a,
                                                       std::function<double(double)>& b, double h, double x){
    return 0.5/(h*h) - 1.5 * a(x) / h;
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FOURTH>::calc(std::function<double(double)>& a,
                                                        std::function<double(double)>& b, double h, double x){
    return -5/(3*(h*h)) + 5/(6*h) * a(x) + b(x);
}

template<>
double Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(std::function<double(double)>& a,
                                                       std::function<double(double)>& b, double h, double x){
    return 11/(12*(h*h)) + 1/(4*h) * a(x);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(std::function<double(double)>& a,
                                                         std::function<double(double)>& b, double h, double x){
    return 1/(h*h) - a(x)/(2*h);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(std::function<double(double)>& a,
                                                          std::function<double(double)>& b, double h, double x){
    return -2/(h*h) + b(x);
}

template<>
double Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(std::function<double(double)>& a,
                                                         std::function<double(double)>& b, double h, double x){
    return 1/(h*h) + a(x)/(2*h);
}

double GetCoef(std::pair<std::string, int> index, std::function<double(double)>& a,
               std::function<double(double)>& b, double h, double x){
    if(index.first == "1")
    {
        if(index.second == 1) return 11/(12*(h*h)) - 0.25*a(x)/h;

        if(index.second == 2) return -5 / (3*(h*h)) - 5/(6*h) * a(x) + b(x);

        if(index.second == 3) return 0.5/(h*h) + 1.5 * a(x) / h;

        if(index.second == 4) return 1/(3*(h*h)) - 0.5/h * a(x);

        if(index.second == 5) return -1/(12*(h*h)) + 1/(12*h) * a(x);
    }

    if(index.first == "i")
    {
        if(index.second == 1) return -1/(12*(h*h)) + a(x)/(12*h);

        if(index.second == 2) return 4/(3 * (h*h)) - 2/(3*h) * a(x);

        if(index.second == 3) return  b(x) - 2.5/(h*h);

        if(index.second == 4) return 4/(3*(h*h)) + 2/(3*h) * a(x);

        if(index.second == 5) return -1/(12*(h*h)) - 1/(12*h) * a(x);
    }

    if(index.first == "n-2")
    {
        if(index.second == 1) return -1/(12*(h*h)) - a(x)*1/(12*h);

        if(index.second == 2) return 1 / (3*(h*h)) + 1/(2*h) * a(x);

        if(index.second == 3) return  0.5/(h*h) - 1.5 * a(x) / h;

        if(index.second == 4) return -5/(3*(h*h)) + 5/(6*h) * a(x) + b(x);

        if(index.second == 5) return 11/(12*(h*h)) + 1/(4*h) * a(x);
    }

    if(index.first == "three diag")
    {
        if(index.second == 1) return 1/(h*h) - a(x)/(2*h);

        if(index.second == 2) return -2/(h*h) + b(x);

        if(index.second == 3) return 1/(h*h) + a(x)/(2*h);
    }

    return 0.;
}

std::pair<Slae::Matrix::FiveDiagonalMatrix, std::vector<double>> ExpandedMatrixForLinearBoundaryValueProblem4(double left_bound_x, double right_bound_x, double left_bound_y,
                                                         double right_bound_y, int number_of_splits,
                                                         std::function<double(double)>& a, std::function<double(double)>& b,
                                                         std::function<double(double)>& f) {
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    Slae::Matrix::FiveDiagonalMatrix data = Slae::Matrix::FiveDiagonalMatrix(number_of_splits + 1);
    data(0, 2) = 1;
    data(data.rows() - 1, 2) = 1;

    //Использование коэффициентов для 3х диагональной матрицы

    data.fill_row(1, 0,
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * 1),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * 1),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * 1), 0);

    data.fill_row(data.rows() - 2, 0,
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * (data.rows() - 2)),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * (data.rows() - 2)),
                  Calc<ROWINDEX::R_THREE_DIAG, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * (data.rows() - 2)), 0);

    //Элементарные преобразования и заполнение коэффициентов 5ти диагональной матрицы

    // 1 < i < n - 2
    for(int i = 2; i < data.rows() - 2; ++i){
        data.fill_row(i, Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * i),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * i),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * i),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h * i),
                      Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * i));
    }
    //i = 1
//    data.fill_row(1, 0, Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h),
//                  Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h));
//
//    data.multiply_row_by_value(1, Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 2) /
//            Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h));
//    data.fill_row(1, 0, data(1, 1) - data(2, 0),
//                  data(1, 2) - data(2, 1),
//                  data(1, 3) - data(2, 2),
//                  data(1, 4) - data(2, 3));
//
//    // i = n - 2
//    data.fill_row(data.rows() - 2, Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h *(data.rows() - 2)),
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h *(data.rows() - 2)) ,
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FOURTH>::calc(a, b, h, left_bound_x + h *(data.rows() - 2)),
//                  Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h *(data.rows() - 2)), 0);
//
//
//    data.multiply_row_by_value(data.rows() - 2,
//                               Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h *(data.rows() - 3)) /
//                                       Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h *(data.rows() - 2)));
//
//    data.fill_row(data.rows() - 2, data(data.rows() - 2, 0) - data(data.rows() - 3, 1),
//                  data(data.rows() - 2, 1) - data(data.rows() - 3, 2),
//                  data(data.rows() - 2, 2) - data(data.rows() - 3, 3),
//                  data(data.rows() - 2, 3) - data(data.rows() - 3, 4), 0);

    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = f(left_bound_x + h * i);
    }

    y.back() = right_bound_y;

    y[1] *= Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 2) /
            Calc<ROWINDEX::R_FIRST, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * 1);

    y[1] -= y[2];

    y[data.rows() - 2] *= Calc<ROWINDEX::R_I, COLUMNINDEX::C_FIRST>::calc( a, b, h,left_bound_x + h * (data.rows() - 3) /
    Calc<ROWINDEX::R_N_2, COLUMNINDEX::C_FIFTH>::calc(a, b, h, left_bound_x + h * (data.rows() - 2)));

    y[data.rows() - 2] -= y[data.rows() - 3];

    return {data, y};
}

#endif //MY_PROJECT_LINEARBOUNDARYVALUEPROBLEM4_HPP
