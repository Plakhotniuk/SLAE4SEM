//
// Created by Арсений Плахотнюк on 23.03.2022.
//

#include "FiveDiagonalMatrix.hpp"
#include <iostream>
using std::cout;
using std::array;
using std::vector;




namespace Slae::Matrix{
    FiveDiagonalMatrix::FiveDiagonalMatrix(unsigned int size): data_(size) {};

    FiveDiagonalMatrix FiveDiagonalMatrix::Zero(unsigned int size)
    {
        FiveDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {0., 0., 0., 0., 0.};
        }
        return result;
    }

    FiveDiagonalMatrix FiveDiagonalMatrix::Identity(unsigned int size)
    {
        FiveDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {0., 0., 1., 0., 0.};
        }
        return result;
    }

    FiveDiagonalMatrix FiveDiagonalMatrix::FiveNumbers(unsigned int size, double a, double b, double c, double d, double e) {
        FiveDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {a, b, c, d, e};
        }
        return result;
    }

    unsigned int FiveDiagonalMatrix::rows() const noexcept
    {
        return data_.size();
    }

    double & FiveDiagonalMatrix::operator()(int i, int j)
    {
#ifndef NDEBUG
        if (i >= data_.size())
        {
            std::stringstream buff;
            buff << "Index i exceeds matrix size! Received index: " << i << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }

        if (j > 4)
        {

            std::stringstream buff;
            buff << "Index j must belong to {0, 1, 2, 3, 4}! Received index: " << i << ". Matrix size: "
                 << data_.size() << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }
#endif
        return data_[i][j];
    }

    const double & FiveDiagonalMatrix::operator()(int i, int j) const
    {
#ifndef NDEBUG
        if (i >= data_.size())
        {
            std::stringstream buff;
            buff << "Index i exceeds matrix size! Received index: " << i << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }

        if (j > 4)
        {
            std::stringstream buff;
            buff << "Index j must belong to {0, 1, 2}! Received index: " << j << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }
#endif
        return data_[i][j];

    }
    void FiveDiagonalMatrix::fill_row(unsigned ind, double a, double b, double c, double d, double e){
        data_[ind][0] = a;
        data_[ind][1] = b;
        data_[ind][2] = c;
        data_[ind][3] = d;
        data_[ind][4] = e;
    }
    DIAG_DOM FiveDiagonalMatrix::check_diagonal_domimance() const {
        bool one_strict_inequality = false;
        bool non_strict_inequality = true;
        for(int i = 0; i < rows(); ++i){
            if(abs(data_[i][2]) <  abs(data_[i][0]) +  abs(data_[i][1]) +
            abs(data_[i][3]) +  abs(data_[i][4])){
                cout << "The sufficient condition of diagonal dominance is not fulfilled in row: " << i << std::endl;
                non_strict_inequality =!non_strict_inequality;
                break;
            }
            if(abs(data_[i][2]) >  abs(data_[i][0]) +  abs(data_[i][1]) +
                                   abs(data_[i][3]) +  abs(data_[i][4])){
                one_strict_inequality = !one_strict_inequality;
            }
        }
        if(one_strict_inequality && non_strict_inequality){
            return DIAG_DOM::NO_DIAG_DOM;
        }
    }
    void FiveDiagonalMatrix::multiply_row_by_value(unsigned int ind, double val) {
        for(int i = 0; i < 5; ++i)
            data_[ind][i] *= val;
    }


}
