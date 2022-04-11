//
// Created by Арсений Плахотнюк on 23.03.2022.
//
#include <iostream>
#include "FiveDiagonalSolver.hpp"
namespace Slae::Solvers {
    std::vector<double>solveFiveDiagonal(const Slae::Matrix::FiveDiagonalMatrix &matrix, const std::vector<double> &col) {
        if (matrix.rows() != col.size()) {
            std::stringstream buff;
            buff << "Matrix and column dimensions are not equal! Matrix size: " << matrix.rows() << ". Column size: "
                 << col.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw SlaeBaseExceptionCpp(buff.str());
        }


        int n = col.size();

        std::vector<double> result(n);
        std::vector<std::array<double, 5>> coefs(n);
        // mu - 0, alfa - 1, betta - 2, z - 3, gamma - 4
        // e - 0; c - 1; d - 2; a - 3; b - 4

        coefs[0][0] = matrix(0, 2);
        coefs[0][1] = matrix(0,3) / coefs[0][0];
        coefs[0][2] = matrix(0, 4) / coefs[0][0];
        coefs[0][3] = col[0] / coefs[0][0];
        coefs[1][4] = matrix(1, 1);
        coefs[1][0] = matrix(1, 2) - coefs[0][1] * coefs[1][4];
        coefs[1][1] = (matrix(1, 3) - coefs[0][2] * coefs[1][4]) / coefs[1][0];
        coefs[1][2] = matrix(1, 4) / coefs[1][0];
        coefs[1][3] = (col[1] - coefs[0][3] * coefs[1][4]) / coefs[1][0];

        for(int i = 2; i < n - 2; ++i){
            coefs[i][4] = matrix(i, 1) - coefs[i - 2][1] * matrix(i, 0);
            coefs[i][0] = matrix(i, 2) - coefs[i - 2][2] * matrix(i, 0) - coefs[i - 1][1] * coefs[i][4];
            coefs[i][1] = (matrix(i, 3) - coefs[i - 1][2] * coefs[i][4]) / coefs[i][0];
            coefs[i][2] = matrix(i, 4) / coefs[i][0];
            coefs[i][3] = (col[i] - coefs[i - 2][3] * matrix(i, 0) - coefs[i - 1][3] * coefs[i][4]) / coefs[i][0];
        }

        coefs[n - 2][4] = matrix(n - 2, 1) - coefs[n - 4][1] * matrix(n - 2, 0);
        coefs[n - 2][0] = matrix(n - 2, 2) - coefs[n-4][2] * matrix(n - 2, 0) - coefs[n - 3][1] * coefs[n - 2][4];
        coefs[n - 2][1] = (matrix(n - 2, 3) - coefs[n-3][2] * coefs[n - 2][4]) / coefs[n - 2][0];
        coefs[n - 1][4] = matrix(n - 1, 1) - coefs[n - 3][1] * matrix(n - 1, 0);
        coefs[n - 1][0] = matrix(n - 1, 2) - coefs[n-3][2] * matrix(n - 1, 0) - coefs[n - 2][1] * coefs[n - 1][4];
        coefs[n - 2][3] = (col[n - 2] - coefs[n - 3][3] * matrix(n - 2, 0) - coefs[n - 3][3] * coefs[n - 2][4]) / coefs[n - 2][0];
        coefs[n - 1][3] = (col[n - 1] - coefs[n - 2][3] * matrix(n - 1, 0) - coefs[n - 2][3] * coefs[n - 1][4]) / coefs[n - 1][0];

        result[n - 1] = coefs[n - 1][3];
        result[n - 2] = coefs[n - 2][3] - coefs[n - 2][1] * result[n - 1];

        for(int i = n - 3; i >= 0; --i){
            result[i] = coefs[i][3] - coefs[i][1] * result[i + 1] - coefs[i][2] * result[i + 2];
        }

        return result;
    }
}