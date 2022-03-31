//
// Created by Арсений Плахотнюк on 23.03.2022.
//

#ifndef MY_PROJECT_FIVEDIAGONALMATRIX_HPP
#define MY_PROJECT_FIVEDIAGONALMATRIX_HPP
#include "my_project/SlaeBaseException.hpp"
#include <vector>
#include <array>
#include <sstream>
#include <string>

namespace Slae::Matrix{
    class FiveDiagonalMatrix{

    private:
        std::vector<std::array<double, 5>> data_;

    public:
        /** @brief ThreeDiagonalMatrix class constructor
         * Creates three-diagonal matrix with size 'size'
         *
         * @param size -- matrix size
        */
        explicit FiveDiagonalMatrix(unsigned int size);
        /** @brief ThreeDiagonalMatrix class static constructor
         * Creates three-diagonal zero matrix with size 'size'
         *
         * @param size -- matrix size
        */
        static FiveDiagonalMatrix Zero(unsigned int size);

        /** @brief ThreeDiagonalMatrix class static constructor
         * Creates three-diagonal identity matrix with size 'size'
         *
         * @param size -- matrix size
        */
        static FiveDiagonalMatrix Identity(unsigned int size);

        /** @brief FiveDiagonalMatrix class static constructor
         * Creates three-diagonal matrix with size 'size' and filled with three numbers
         *
         * @param a -- number of lower diagonal
         * @param b -- number of main diagonal
         * @param c -- number of upper diagonal
         * @param size -- matrix size
        */
        static FiveDiagonalMatrix FiveNumbers(unsigned int size, double a, double b, double c, double d, double e);
        /** @brief FiveDiagonalMatrix class static constructor
         * Creates three-diagonal matrix with size 'size' and filled with three numbers
         *
         * @param i -- index of row
         * @param j -- index of column
         * @param size -- matrix size
        */
        double & operator()(int i, int j);

        /** @brief FiveDiagonalMatrix class method
        *
        * @param i -- index of row
        * @param j -- index of column
        * @return matrix element with i,j indexes
        */
        [[nodiscard]] const double & operator()(int i, int j) const;
        /** @brief FiveDiagonalMatrix class method
         * @return matrix row-size
        */
        [[nodiscard]] unsigned int rows() const noexcept;
        /**
         * @brief FiveDiagonalMatrix class method
         * @param ind -- index of row
         * @param  a -- number of lower diagonal
         * @param b -- number of main diagonal
         * @param c -- number of upper diagonal
         */
        void fill_row(unsigned ind, double a, double b, double c, double d, double e);

        void multiply_row_by_value(unsigned ind, double val);

        void substract_row1_from_row2(unsigned ind1, unsigned ind2);

        void check_diagonal_domimance() const;

    };

}   //namespace Slae::Matrix

#endif //MY_PROJECT_FIVEDIAGONALMATRIX_HPP
