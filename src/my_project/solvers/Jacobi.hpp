//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_JACOBI_HPP
#define SLAE_JACOBI_HPP
#include "../sparse/CSR.hpp"
#include "my_project/utility/Norm.hpp"
template<typename T>
/**
 * функция, реализующая метод Якоби для решения СЛАУ
 * @tparam T - шаблонный тип
 * @param A - матрица системы
 * @param b - свободный столбец системы
 * @param initialState - начальное приближение решения
 * @param tolerance - точность, предъявляемая к результату
 * @return currentState - решение в пределах tolarance
 * @throw SlaeBaseException - выбрасывается при неправильных размерностях входных данных
 */
std::vector<T> Jacobi(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
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
            tempState[i] = (b[i] - sum)/ A(i, i);
        }
        currentState = tempState;
        r = A * currentState - b;
    }
    return currentState;
};

#endif//SLAE_JACOBI_HPP
