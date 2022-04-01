//
// Created by Арсений Плахотнюк on 11.03.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include <my_project/solvers/NonlinearBoundaryValueProblem2.hpp>
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>
#include <algorithm>


TEST(NONLINEARBOUNDARYVALUEPROBLEM, TEST) {
    std::function<double(double)> f = [](double x) { return 2*x; };

    std::function<double(double)> a = [](double x) { return 0.; };

    std::function<double(double)> b = [](double x) { return -1.; };

    std::function<double(double)> y_func = [](double x) { return sinh(x) / sinh(1.) - 2*x; };

    double left_bound_x = 0.;
    double right_bound_x = 1.;
    double left_bound_y = 0.;
    double right_bound_y = -1.;

    int max_number_of_splits = 220;
    std::vector<double> y(max_number_of_splits + 1);
    std::fstream file;
    file.open("test_2_func3.txt", std::fstream::out);
       std::vector<double> x(max_number_of_splits + 1);
        double h = (right_bound_x - left_bound_x)/max_number_of_splits;
        for(int i = 0; i < x.size(); ++i){
            x[i] = left_bound_x + h * i;
            y[i] = y_func(x[i]);
        }
    std::vector<double> solution = NonlinearBoundaryValueProblem2(left_bound_x, right_bound_x,
                                                                      left_bound_y, right_bound_y,
                                                                     max_number_of_splits, a, b, f);
    for(int i=0; i < solution.size(); ++i){
        file << solution[i] << " " << x[i] << " "<<y[i];
        file << '\n';
    }
    file.close();
}

TEST(NONLINEARBOUNDARYVALUEPROBLEM2, TEST2) {
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
    file.open("test_2_err3.txt", std::fstream::out);
    for(int j = 20; j < max_number_of_splits; j += 10){
        std::vector<double> solution = NonlinearBoundaryValueProblem2(left_bound_x, right_bound_x,
                                                                      left_bound_y, right_bound_y,
                                                                      j, a, b, f);
        std::vector<double> x(j + 1);
        std::vector<double> y(j + 1);

        double h = (right_bound_x - left_bound_x)/j;
        for(int i = 0; i < x.size(); ++i){
            x[i] = left_bound_x + h * i;
            y[i] = y_func(x[i]);
        }
        std::vector<double> err(j);
        for(int i=0; i<j; i++){
            err[i] = abs(solution[i] - y[i]);
        }
        double max = *std::max_element(err.begin(), err.end());
        file << j << " "<< max << " ";
        file << '\n';
    }
    file.close();

}



