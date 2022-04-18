//
// Created by Арсений Плахотнюк on 19.03.2022.
//

#ifndef MY_PROJECT_STEEPEST_GRADIENT_DESCENT_HPP
#define MY_PROJECT_STEEPEST_GRADIENT_DESCENT_HPP


#include "../sparse/CSR.hpp"
#include "../utility/Norm.hpp"
#include <my_project/utility/Overloads.hpp>
#include <sstream>
#include "my_project/Exceptions/SlaeBaseException.hpp"

template <typename T>
std::vector<T> Steepest_gradient_descent(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &init, const T &tolerance) {
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
    if (init.size() != b.size()) {
        std::stringstream buff;
        buff << "InitState and free column must have be same size! init.size " << init.size() << "! b.size << "
             << b.size() << "!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
    std::vector<T> res = init;
    std::vector<T> r = A * res - b;
    T tao = (r * r) / (r * (A * r));
    int i = 0;
    while (norm(r, NormType::ThirdNorm) > tolerance) {
        res = res - tao * r;
        r = A * res - b;
        tao = (r * r) / (r * (A * r));
        ++i;
    }
    return res;
#endif // NDEBUG
}

#endif //MY_PROJECT_STEEPEST_GRADIENT_DESCENT_HPP
