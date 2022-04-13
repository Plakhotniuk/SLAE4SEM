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
std::vector<T> Gmres(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
    std::vector<T> r = A * initialState - b;
    std::vector<T> currentState = initialState;
    std::vector<std::vector<T>> v(currentState.size() + 1);
    DenseMatrix<T> h(currentState.size() + 1, currentState.size());
    DenseMatrix<T> R(currentState.size() + 1, currentState.size());
    std::vector<std::vector<T>> y(currentState.size() + 1);
    int N;
    int m = 6;
    bool solution = true;
    std::vector<std::vector<T>> z(currentState.size() + 1);
    double gamma;
    int j = 0;
    while(solution){
        std::cout<<"j: "<<j<<std::endl;
        double beta = norm(r, NormType::SecondNorm);
        v[0] = r / beta;
        std::vector<double> e(currentState.size() + 1);
        e[0] = 1.;
        z[0] = beta * e;
        for(int i = 0; i < m; ++i)
        {
            std::cout<<"i: "<<i<<std::endl;
            if(solution)
            {
                std::vector<T> t = A * v[i];
                for(int k = 0; k < i; ++k)
                {
                    h(k,i) = t * v[k];
                    t = t - h(k,i) * v[k];
                }
                h(i+1, i) = norm(t, NormType::SecondNorm);
                v[i+1] = t / h(i+1, i);
                R(0, i) = h(0, i);
                std::vector<double> c(m);
                std::vector<double> s(m);
                for(int k = 1; k < i; ++k)
                {
                    std::cout<<"k: "<<i<<std::endl;
                    double alpha = c[k-1] * R(k-1, i) + s[k-1] * h(k, i);
                    R(k, i) = c[k-1] * h(k, i) - s[k-1] * R(k-1, i);
                    R(k-1, i) = alpha;
                }
                c[i] = R(i, i) / sqrt(R(i, i)*R(i, i) + h(i+1, i)*h(i+1, i));;
                s[i] = h(i+1, i) / sqrt(R(i, i)*R(i, i) + h(i+1, i)*h(i+1, i));;
                z[i+1] = -s[i] * z[i];
                z[i] = c[i] * z[i];
                gamma = norm(z[i+1], NormType::SecondNorm);
                if(gamma < tolerance)
                {
                    N = i;
                    solution = !solution;
                }
            }
        }
        if(solution)
        {
            N = m;
            y[N] = z[N] / R(N, N);
            solution = !solution;
        }
        if(!solution)
        {
            for(int k = N - 1; k > 0; k--)
            {
                std::vector<T> sum1(y.size());
                for(int i = k + 1; i < N; ++i)
                {
                    sum1 = sum1 + R(k, i) * y[i];
                }
                y[k] = z[k] - sum1 / R(k, k);

                std::vector<T> sum2(v[0].size());
                for(int i = 0 ; i < N; ++i)
                {
                    sum2 = sum2 + v[i] * y[i];
                }
                currentState = currentState - sum2;
            }
            if(abs(gamma) < tolerance)
                return currentState;
            r = A * currentState - b;
            solution = !solution;
        }
        ++j;
    }
}