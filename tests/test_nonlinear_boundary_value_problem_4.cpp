//
// Created by Арсений Плахотнюк on 24.03.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include <my_project/solvers/NonlinearBoundaryValueProblem4.hpp>
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>
#include <algorithm>


TEST(NONLINEARBOUNDARYVALUEPROBLEM2, TEST) {
    std::function<double(double)> f = [](double x) { return 1.; };

    std::function<double(double)> a = [](double x) { return 0.; };

    std::function<double(double)> b = [](double x) { return 1.; };

    std::function<double(double)> y_func = [](double x) { return cos(x); };

    double left_bound_x = 0.;
    double right_bound_x = M_PI;
    double left_bound_y = 1.;
    double right_bound_y = -1.;

    int max_number_of_splits = 20;
    std::fstream file;
    file.open("test_4.txt", std::fstream::out);

    std::vector<double> solution = NonlinearBoundaryValueProblem4(left_bound_x, right_bound_x,
                                                                  left_bound_y, right_bound_y,
                                                                  max_number_of_splits, a, b, f);
    std::vector<double> x(max_number_of_splits + 1);
    std::vector<double> y(max_number_of_splits + 1);

    double h = (right_bound_x - left_bound_x)/max_number_of_splits;
    for(int i = 0; i < x.size(); ++i){
        x[i] = left_bound_x + h * i;
        y[i] = y_func(x[i]);
    }
    std::vector<double> err(max_number_of_splits);
    std::cout << "[";
    for(int i=0; i<max_number_of_splits; i++){
        err[i] = abs(solution[i] - y[i]);
        std::cout << solution[i] << " ";
    }
    std::cout << "]";
    std::cout<< std::endl;
    double max = *std::max_element(err.begin(), err.end());
    file << max_number_of_splits << " "<< max << " ";
    file << '\n';

    
    file.close();
}


