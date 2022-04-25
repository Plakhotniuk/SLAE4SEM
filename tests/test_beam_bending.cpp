//
// Created by Арсений Плахотнюк on 18.04.2022.
//
#include "gtest/gtest.h"
#include <my_project/utility/Overloads.hpp>
#include "my_project/matrix/NonLinearBoundaryValueProblemMatrix5.hpp"
#include <iostream>
#include "functional"
#include <cmath>
#include <fstream>

TEST(BEAMBENDING, BEAM) {
    int number_of_loops = 200;
    std::fstream file;
    file.open("test_beam6.txt", std::fstream::out);

    for(int q = 0; q < number_of_loops; q++) {

        std::function<double(double, double)> f = [q](double x, double y) {
            return -(0.03*q)*(0.03*q) * sin(y);
        };

        std::function<double(double, double, double)> a = [](double x, double y, double y_dev) { return 0.; };

        std::function<double(double, double, double)> b = [](double x, double y, double y_dev) { return 0.; };

        std::function<double(double)> y_func = [](double x) { return asin(x); };

        double left_bound_x = 0.;
        double right_bound_x = 1;
        double left_bound_y = 0.;
        double right_bound_y = M_PI;

        int max_number_of_splits = 300;
        int n_iterations = 30;

        double h = (right_bound_x - left_bound_x) / max_number_of_splits;

        std::pair<std::vector<double>, std::vector<double>> solution;

        std::pair<std::vector<double>, std::vector<double>> y_0_y_0_der = InitialApproach_y5(left_bound_x,
                                                                                             right_bound_x,
                                                                                             left_bound_y,
                                                                                             right_bound_y,
                                                                                             max_number_of_splits,
                                                                                             y_func);

        ///approach
        solution = y_0_y_0_der;
        for (int i = 1; i < n_iterations; ++i) {
            std::pair<Slae::Matrix::FiveDiagonalMatrix, std::vector<double>> matrix =
                    ExpMatrixNonLinearBVP5(left_bound_x, right_bound_x,
                                           left_bound_y, right_bound_y,
                                           max_number_of_splits, a, b, f, solution);
            solution.first = Slae::Solvers::solveFiveDiagonal(matrix.first, matrix.second);
        }

        ///Write x to file
        std::vector<double> x(solution.first.size());
        file << x[0] << " ";
        for (int i = 1; i < solution.first.size(); ++i) {
            x[i] = x[i - 1] + h * cos(solution.first[i]);
            file << x[i] << " ";
        }

//        for (int i = 0; i < solution.first.size() - 1; ++i) {
//            x.push_back(2*x[solution.first.size() - 1] - x[solution.first.size() - 2 - i]);
//            file << x[i + solution.first.size()] << " ";
//        }
        file << "\n";

        ///Write y to file
        std::vector<double> y_x(solution.first.size());
        file << y_x[0] << " ";
        for (int i = 1; i < solution.first.size(); ++i) {
            y_x[i] = y_x[i - 1] + h * sin(solution.first[i]);
            file << y_x[i] << " ";
        }

//        for (int i = 0; i < solution.first.size() - 1; ++i) {
//            y_x.push_back(y_x[solution.first.size() - 2 - i]);
//            file << y_x[i + solution.first.size()] << " ";
//        }

        file << "\n";

        ///Write l and M(x)/(EI) to file
        for (int i = 0; i < y_x.size(); ++i) {
            file << left_bound_x + h * i << " ";
        }
        file << "\n";

        std::vector<double> y_second_der = Approach_y_second_der5(y_x, h);
        std::vector<double> y_der = Approach_y_der5(y_x, h);
        for (int i = 0; i < y_x.size(); ++i) {
            file << y_second_der[i] / pow(1 + y_der*y_der, 1.5) << " ";
        }
        file << "\n";

    }
    file.close();

}

