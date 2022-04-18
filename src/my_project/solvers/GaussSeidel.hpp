//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_GAUSSSEIDEL_HPP
#define SLAE_GAUSSSEIDEL_HPP

#include <my_project/utility/Norm.hpp>
#include "../sparse/CSR.hpp"
#include <sstream>
#include "my_project/Exceptions/SlaeBaseException.hpp"



template<typename T>
std::vector<T> GaussSeidel(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
#ifndef NDEBUG
    if (A.sizeH() != A.sizeW()) {
        std::stringstream buff;
        buff << "Matrix must be quadratic!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
    if (A.sizeH() != b.size()) {
        std::stringstream buff;
        buff << "Matrix and free column must have be same size! A.size " << A.sizeH() << "! b.size << " << b.size()
             << "!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
    if (initialState.size() != b.size()) {
        std::stringstream buff;
        buff << "InitState and free column must have be same size! initialState.size " << initialState.size() << "! b.size << "
             << b.size() << "!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
#endif // NDEBUG
    std::vector<T> r = A * initialState - b;
    std::vector<T> currentState = initialState;
    std::vector<T> tempState(currentState.size());
    T sum;
    while(norm(r, NormType::ThirdNorm) > tolerance){
        for(int i = 0; i < A.sizeH(); ++i){
            sum = static_cast<T>(0);
            int skip = A.rows[i];
            int count = A.rows[i+1] - skip;
            for(int k = skip; k < skip + count; ++k){
                if(A.cols[k] != i) sum += A.values[k] * currentState[i];
            }
            currentState[i] = (b[i] - sum)/ A(i, i);
        }
        r = A * currentState - b;
    }
    return currentState;
};

#endif//SLAE_GAUSSSEIDEL_HPP
