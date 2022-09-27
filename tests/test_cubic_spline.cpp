//
// Created by Арсений Плахотнюк on 27.09.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include "my_project/solvers/ThreeDiadonalSolver.hpp"
#include <iostream>
#include "functional"
#include <fstream>

double get_dif(const std::vector<double> &x, const std::vector<double> &f)
{
    if (f.size() > 2)
    {
        std::vector<double> f_left(f.size() - 1);
        std::vector<double> x_left(x.size() - 1);
        std::vector<double> f_right(f.size() - 1);
        std::vector<double> x_right(x.size() - 1);
        for(int i = 0; i < f.size() - 1; ++i){
            f_left[i] = f[i];
            x_left[i] = x[i];
            f_right[i] = f[i+1];
            x_right[i] = x[i+1];
        }
        return (get_dif(x_left, f_left) - get_dif(x_right, f_right)) / (x[x.size()-1] - x[0]);
    }
    else if (f.size() == 2) {
        return (f[1] - f[0]) / (x[1] - x[0]);
    }
}

TEST(THREEDIAGSPLINE, TEST) {

    std::vector<double> x = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> y = {0.0, 0.033, 0.067, 0.1, 0.134, 0.168, 0.203, 0.238, 0.273, 0.309, 0.346};
    int n = x.size();
    std::vector<double> h(n-1);

    for(int i = 0; i < n-1; ++i){
        h[i] = x[i + 1] - x[i];
    }

    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(n - 2);
    std::vector<double> u(n-2);
    for(int i = 1; i < n-1; ++i){
        u[i-1] = 6 * get_dif(std::vector({x[i-1], x[i], x[i+1]}), std::vector({y[i-1], y[i], y[i+1]}));
    }
    matrix.fill_row(0, 0, 2, h[1] / (h[0] + h[1]));

    for(int i = 1; i < n - 3; ++i){
        matrix.fill_row(i,
                        h[i-1] / (h[i-1] + h[i]),
                        2,
                        h[i-1] / (h[i-1] + h[i]));
    }
    matrix.fill_row(n - 3, h[n-3] / (h[n-3] + h[n-2]), 2, 0);

    std::vector<double> c =  Slae::Solvers::solveThreeDiagonal(matrix, u);
    c.insert(c.begin(), 0);
    c.emplace_back(0);

    std::vector<double> a = y;
    std::vector<double> b(n);

    for (int i = 0; i < n; ++i){
        b[i] = (c[i+1] * h[i+1] / 3) + (c[i] * h[i+1] / 6) + get_dif(x, {y[i], y[i+1]});
    }

    std::vector<double> d(n);
    d[0] = c[1] / h[0];
    for (int i = 1; i < n; ++i){
        d[i] = (c[i+1] - c[i]) / h[i];
    }
    double var = 0.95;
    int i = 10;
    for(int j = 0; j < n - 1; ++j){
        if(x[j] <= var && var < x[j+1]){
            i = j;
        }
    }
    double res = a[i] + b[i] * (var - x[i]) +
            c[i] * (var - x[i])*(var - x[i]) / 2 + d[i] *(var - x[i])*(var - x[i])*(var - x[i])/6;
    std::cout <<std::setprecision(9) << res <<std::endl;
}




