//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_NORM_HPP
#define SLAE_NORM_HPP
#include "vector"
#include "cmath"

enum class NormType{
    FirstNorm = 1,
    SecondNorm = 2,
    ThirdNorm = 3,
};
/**
 * функция, реализующая различные нормы вектора
 * @tparam T шаблонный тип
 * @param vector вектор
 * @param normType тип нормы
 * @return norm соответствующая норма вектора
 */
template <typename T>
const T& norm(const std::vector<T>& vector, NormType normType){
    T norm = static_cast<T>(0);
    if(normType == NormType::FirstNorm){
        for(const auto& elm: vector){
            T abs = std::abs(elm);
            if(abs > norm) norm = abs;
        }
    }
    if(normType == NormType::SecondNorm){
        for(const auto& elm: vector){
            norm += std::abs(elm);
        }
    }
    if(normType == NormType::ThirdNorm){
        for(const auto& elm: vector){
            norm += elm*elm;
        }
        norm = std::sqrt(norm);
    }
    return norm;
}


#endif//SLAE_NORM_HPP
