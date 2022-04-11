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
    DenseMatrix<T> h(currentState.size(), currentState.size() + 1);
    std::vector<std::vector<T>> R(currentState.size());
    std::vector<T> y(currentState.size());
    int N;
    int m = 4;
    bool solution = true;
    std::vector<T> z;
    double gamma;
    int j = 0;
    while(solution){
        double beta = norm(r, NormType::SecondNorm);
        v[0] = r / beta;
        std::vector<double> e;
        e[0] = 1.;
        z = beta * e;
        for(int i = 0; i < m; ++i)
        {
            if(solution)
            {
                std::vector<T> t = A * v[i];
                for(int k = 0; k < i; ++k)
                {
                    h(k,j) = t * v[k];
                    t = t - h(k,j) * v[k];
                }
                h(i+1, i) = norm(t, NormType::SecondNorm);
                v[i+1] = t / h(i+1, i);
                R[0][i] = h(0, i);
                double delta;
                std::vector<double> c(m);
                std::vector<double> s(m);
                for(int k = 1; k < i; ++k)
                {
                    double alpha = c[k-1] * R[k-1][i] + s[k-1] * h(k, i);
                    R[k][i] = c[k-1] * h(k, i) - s[k-1] * R[k-1][i];
                    R[k-1][i] = alpha;
                }
                delta = sqrt(R[i][i]*R[i][i] + h(i+1, i)*h(i+1, i));
                c[i] = R[i][i] / delta;
                s[i] = h(i+1, i) / delta;
                z[i+1] = -s[i] * z[i];
                z[i] = c[i] * z[i];
                gamma = abs(z[i+1]);
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
            y[N] = z[N] / R[N][N];
            solution = !solution;
        }
        if(!solution)
        {
            for(int k = N - 1; k > 0; k--)
            {
                double sum1 = 0.;
                for(int i = k + 1; i < N; ++i)
                {
                    sum1 = sum1 + R[k][i] * y[i];
                }
                y[k] = z[k] - sum1 / R[k][k];
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