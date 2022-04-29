//
// Created by Арсений Плахотнюк on 09.04.2022.
//

#ifndef MY_PROJECT_GMRES_H
#define MY_PROJECT_GMRES_H


#include "../sparse/CSR.hpp"
#include "../dense/Densematrix.hpp"
#include "my_project/utility/Norm.hpp"
#include "my_project/utility/Overloads.hpp"
#include <ostream>

template<typename T>
void AddVecToKrylovSubspace(const CSR<T>& A, DenseMatrix<T>& V, DenseMatrix<T>& H, int j){
    std::vector<T> temp(H.sizeW());
    temp = A * V.get_col(j);

    for(int k = 0; k < j + 1; ++k)
    {
        H(k, j) = V.get_col(k) * temp;
        temp = temp - H(k, j)*V.get_col(k);
    }

    H(j + 1, j) = Norm<T, NormType::ThirdNorm>::get_norm(temp);
    V.write_col(temp / H(j + 1, j), j + 1);
}


template<typename T>
void GivensRotation(DenseMatrix<T>& H, std::vector<T>& z, std::vector<std::pair<T, T>>& rotate_cos_sin, int k){
    double cos;
    double sin;
    double alpha;
    double beta;
    for (int i = 0; i < k; i++)
    {
        beta = H(i, k);
        H(i, k) = H(i, k) * rotate_cos_sin[i].first - H(i + 1, k) * rotate_cos_sin[i].second;
        H(i + 1, k) = beta * rotate_cos_sin[i].second + H(i + 1, k) * rotate_cos_sin[i].first;
    }

    alpha = std::sqrt(H(k, k) * H(k, k) + H(k + 1, k) * H(k + 1, k));
    cos = H(k, k) / alpha;
    sin = - H(k + 1, k) / alpha;

    std::pair<double, double> a(cos, sin);
    rotate_cos_sin[k] = a;
    beta = H(k, k);
    H(k, k) = H(k, k) * cos - H(k + 1, k) * sin;
    H(k + 1, k) = beta * sin + H(k + 1, k) * cos;

    beta = z[k];
    z[k] = z[k] * cos - z[k + 1] * sin;
    z[k + 1] = beta * sin + z[k + 1] * cos;
}

template <typename T>
std::vector<T> ReverseGauss(const DenseMatrix<T>& R, const std::vector<T>& b, const int& n)
{
    std::vector<T> sol(n, 0);

    for (int i = 0; i < n; i++)
    {
        sol[n - i - 1] = b[n - i - 1];
        for (int k = n - 1; k > n - i - 1; k--)
            sol[n - i - 1] -=  R(n - i - 1, k) * sol[k];
        if (R(n - i - 1, n - i - 1) != 0.)
            sol[n - i - 1] /= R(n - i - 1, n - i - 1);
    }

    return sol;
}

template<typename T>
std::vector<T> GMRES(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
    int i = A.sizeH();
    std::vector<T> r = A * initialState - b;
    std::vector<T> currentState = initialState;
    DenseMatrix<T> V(i, i + 1); // ОНБ в подпространстве крылова
    DenseMatrix<T> H(i + 1, i); // матрица Хессенберга
    bool finished = false;
    std::vector<T> z(i + 1);
    std::vector<std::pair<T, T>> rotate_cos_sin(i + 1, {0., 0.});

    while (!finished)
    {
        std::fill(z.begin(), z.end(), 0.);
        z[0] = Norm<T, ThirdNorm>::get_norm(r);
        V.write_col(r, 0);
        if (Norm<T, ThirdNorm>::get_norm(r) != 0.)
            V.write_col(r / Norm<T, ThirdNorm>::get_norm(r), 0);

        for(int j = 1; j < i + 1; ++j)
        {
            AddVecToKrylovSubspace(A, V, H, j - 1);
            GivensRotation(H, z, rotate_cos_sin, j - 1);

            if (std::abs(z[j]) < tolerance)
            {
                std::vector<T> y = ReverseGauss(H, z, j);
                for (int k = 0; k < j; k++)
                    currentState = currentState -  y[k] * V.get_col(k);
                finished = true;
                break;
            }
        }
        if (!finished)
        {
            std::vector<double> y = ReverseGauss(H, z, i);
            for (int k = 0; k < i; k++)
                currentState = currentState -  y[k] * V.get_col(k);
        }
        r = A * currentState - b;

    }
    return currentState;
}

#endif //MY_PROJECT_GMRES_H
