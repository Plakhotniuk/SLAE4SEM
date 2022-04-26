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
    V.write_col(A * V.get_col(j), j + 1);

    for(int k = 0; k < j + 1; ++k)
    {
        H(k, j) = V.get_col(k) * V.get_col(j + 1);
        V.write_col(V.get_col(j + 1) - H(k,j) * V.get_col(k), j + 1);
    }

    H(j + 1, j) = norm(V.get_col(j + 1), NormType::SecondNorm);
    V.write_col(V.get_col(j) / H(j + 1, j), j + 1);
}


template<typename T>
void GivensRotation(DenseMatrix<T>& H, std::vector<T>& z, std::vector<std::pair<T, T>>& rotate_cos_sin, int j){
    rotate_cos_sin[j].first = H(j, j) / sqrt(H(j, j)* H(j, j) + H(j + 1, j)* H(j + 1, j));
    rotate_cos_sin[j].second = -H(j + 1, j) / sqrt(H(j, j)* H(j, j) + H(j + 1, j)* H(j + 1, j));

    for(int i = 0; i < j + 1; ++i)
    {
        H(i, j) = H(i, j) * rotate_cos_sin[i].first - H(i + 1, j) * rotate_cos_sin[i].second;
        H(i + 1, j) = H(i, j) * rotate_cos_sin[i].second +  H(i + 1, j) * rotate_cos_sin[i].first;
    }
    z[j] = z[j] * rotate_cos_sin[j].first -  z[j + 1] * rotate_cos_sin[j].second;
    z[j + 1] = z[j] * rotate_cos_sin[j].second + z[j + 1] * rotate_cos_sin[j].first;
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
    std::vector<T> z(i + 1, 0);
    std::vector<std::pair<T, T>> rotate_cos_sin(i + 1);
    int n = 0;

    while (!finished)
    {
        z[0] = norm(r, NormType::SecondNorm);
        V.write_col(r / norm(r, NormType::SecondNorm), 0);
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
            std::cout<<n<<std::endl;
            ++n;
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
