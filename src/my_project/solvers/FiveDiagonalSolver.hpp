//
// Created by Арсений Плахотнюк on 23.03.2022.
//

#ifndef MY_PROJECT_FIVEDIAGONALSOLVER_HPP
#define MY_PROJECT_FIVEDIAGONALSOLVER_HPP
#include "my_project/matrix/FiveDiagonalMatrix.hpp"
#include "my_project/SlaeBaseException.hpp"
#include <sstream>
#include <vector>

namespace Slae::Solvers
{
    /** @brief Метод решает систему уравнений при помощи метода прогонки
* Решает систему линейных алгебраических уравнений при помощи метода прогонки. О методе прогонки можно узнать из
* https://mipt.ru/education/chair/computational_mathematics/study/materials/compmath/other/Aristova_Zavyalova_Lobanov_2014.pdf
*
* @param matrix трехдиагональная матрица
* @param col столбец правой части
* @return решение СЛАУ
*
* @throw SlaeBaseExceptionCpp выбрасывается, если количество строк матрицы и высота столбца не совпадают
*/
    [[nodiscard]] std::vector<double> solveFiveDiagonal(const Matrix::FiveDiagonalMatrix &matrix, const std::vector<double> &col);

}  // namespace Slae::Solvers
#endif //MY_PROJECT_FIVEDIAGONALSOLVER_HPP
