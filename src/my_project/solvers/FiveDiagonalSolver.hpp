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
* @param matrix пятидиагональная матрица
* @param col столбец правой части
* @return решение СЛАУ
*
* @throw SlaeBaseExceptionCpp выбрасывается, если количество строк матрицы и высота столбца не совпадают
*/
    [[nodiscard]] std::vector<double> solveFiveDiagonal(const Matrix::FiveDiagonalMatrix &matrix, const std::vector<double> &col);

}  // namespace Slae::Solvers
#endif //MY_PROJECT_FIVEDIAGONALSOLVER_HPP
