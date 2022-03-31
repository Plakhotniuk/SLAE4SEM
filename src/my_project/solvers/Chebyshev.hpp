//
// Created by Арсений Плахотнюк on 05.03.2022.
//
#include <cmath>
#include <vector>
#ifndef MY_PROJECT_CHEBYSHEV_H
#define MY_PROJECT_CHEBYSHEV_H
template<typename T, int powOf2>
class ChebyshevPolynomial{
public:
    /**
     * Метод нахождения корней полинома Чебышева
     * @tparam T - шаблонный параметр
     * @param section
     * @return result - вектор решений
     */
    std::vector<T> CalculateRoots(std::pair<T,T> section){
        int polyOrder = std::pow(2, powOf2);
        std::vector<T> roots(polyOrder);
        for(int i = 1; i <= polyOrder; ++i){
            roots[i-1] = (section.first + section.second) / 2. +
                    (section.second - section.first) / 2. +
                    static_cast<T>(std::cos((2 * i - 1) * M_PI_2 / polyOrder));
        }
        std::vector<int> idx = {0, 1};
        std::vector<int> next_ind;
        int currOrder;
        for(int i = 1; i < powOf2; ++i){
            currOrder = std::pow(2, i + 1);
            next_ind.resize(currOrder);
            for(int j = 0; j < currOrder - 1; j += 2){
                next_ind[j] = idx[j/2];
                next_ind[j+1] = currOrder - 1 - next_ind[j];
            }
            idx = next_ind;
        }
        std::vector<T> result(roots.size());
        for(int i = 0; i < result.size(); ++i){
            result[i] = roots[idx[i]];
        }
        return result;
    }
};

#endif //MY_PROJECT_CHEBYSHEV_H
