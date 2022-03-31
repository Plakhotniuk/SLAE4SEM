//
// Created by Арсений Плахотнюк on 09.02.2022.
//
#include "gtest/gtest.h"
#include "my_project/solvers/FiveDiagonalSolver.hpp"
#include "my_project/matrix/FiveDiagonalMatrix.hpp"
#include<iostream>
TEST(TEST_FIVEDIAGMATRIX, NODISCARD){
    auto a = Slae::Matrix::FiveDiagonalMatrix::FiveNumbers(5, 1., 1., 10., 1., 1.);
    std::vector<double> b = {10, 10, 10, 10, 10};
    std::vector<double> solution = Slae::Solvers::solveFiveDiagonal(a, b);
    std::vector<double> online_solver = {0.85551331, 0.769961, 0.674904, 0.769961, 0.85551331};
    for(int i = 0; i < solution.size(); ++i){
        std::cout << solution[i] << " " <<  online_solver[i] << std::endl;
        ASSERT_NEAR(solution[i], online_solver[i], 0.01);
    }
}
