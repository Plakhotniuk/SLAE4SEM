//
// Created by Арсений Плахотнюк on 26.03.2022.
//

#ifndef MY_PROJECT_CG_HPP
#define MY_PROJECT_CG_HPP

#include <my_project/utility/Norm.hpp>
#include "../sparse/CSR.hpp"
#include <sstream>
#include <my_project/SlaeBaseException.hpp>
#include <fstream>
#include <algorithm>

template<typename T>
std::vector<T> CG(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState, const T &tolerance){
    std::vector<T> r = A * initialState - b;
    std::vector<T> currentState = initialState;
    std::vector<T> d = r;
    std::vector<T> tempState(currentState.size());
//    int i = 0;
    T sum;
    std::fstream file;
    file.open("cg_steps.txt", std::fstream::out);
    while(norm(r, NormType::ThirdNorm) > tolerance){
        for(int j = 0; j < currentState.size(); ++j)
            file << currentState[j] << " ";
        file << '\n';
        currentState = currentState - (r * r)/(d * (A * d)) * d;
        tempState = r;
        r = A * currentState - b;
        if(norm(r, NormType::ThirdNorm) < tolerance)
            break;
        else
            d = r + (r * r) / (tempState * tempState) * d;
//        ++i;
    }
    for(int j = 0; j < currentState.size(); ++j)
        file << currentState[j] << " ";
    file << '\n';
    file.close();


    return currentState;
};
#endif //MY_PROJECT_CG_HPP
