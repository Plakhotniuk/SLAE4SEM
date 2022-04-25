//
// Created by Арсений Плахотнюк on 09.04.2022.
//

#ifndef MY_PROJECT_GMRES_H
#define MY_PROJECT_GMRES_H

#endif //MY_PROJECT_GMRES_H
#include "../sparse/CSR.hpp"
#include "../dense/Densematrix.hpp"
#include "my_project/utility/Norm.hpp"
#include "my_project/utility/Overloads.hpp"
#include <ostream>

template<typename T>
std::vector<T> ArnoldiOrthogonalization(std::vector<T>& r_0, const CSR<T> &A){
    std::vector<T> V(A.sizeH());
    std::vector<T> temp(A.sizeH());
    DenseMatrix<T> H(A.size() + 1, A.size());
    V[0] = r_0 / norm(r_0, NormType::SecondNorm);
    for(int i = 1; i < A.sizeH(); ++i)
    {
        temp = A * V[i];
        for(int k = 0; k < i ; ++k)
        {
            H[k][i] = V[k] * temp;
            temp = temp - H[k][i] * V[k];
        }
        H[i+1][i] = norm(temp, NormType::SecondNorm);
        V[i] = temp / H[i+1][i];
    }
    return V;
}

template<typename T>
std::vector<T> GivensRotation(DenseMatrix<T> H, std::vector<T> b){
    for(int i = 0; i < H.sizeW(); ++i)
    {

    }
}

template<typename T>
std::vector<T> GMRES(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
    std::vector<T> r = A * initialState - b;
    std::vector<T> currentState = initialState;
    std::vector<std::vector<T>> V(currentState.size()); // ОНБ в подпространстве крылова
    std::vector<std::vector<T>> y(currentState.size()); // координаты в базисе V
    DenseMatrix<T> h(currentState.size() + 1, currentState.size());
    DenseMatrix<T> R(currentState.size() + 1, currentState.size());
    int m = 4;
    for(int i = 0; i < m; ++i)
    {

    }
}