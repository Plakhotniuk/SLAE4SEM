//
// Created by kuznetsov on 04.02.2022.
//

#include "gtest/gtest.h"
#include <my_project/SlaeBaseException.hpp>
#include "my_project/solvers/ThreeDiadonalSolver.hpp"

void throwException() { throw Slae::SlaeBaseExceptionCpp("Hi"); }

TEST(THREEDIAGONALMATRIX, TEST1) {
    int n = 3;
    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::ThreeNumbers(n, 2.0, 3.0, 1.0);

    std::vector<double> col = {1., 2., 3.};

    std::vector<double> solution = Slae::Solvers::solveThreeDiagonal(matrix, col);

    //got solution using online solver
    std::vector<double> online_solution = {4./15, 0.2, 13./15.};

    for(int i = 0; i < n; ++i){
        ASSERT_NEAR(solution[i], online_solution[i], 0.0001);
    }

}
