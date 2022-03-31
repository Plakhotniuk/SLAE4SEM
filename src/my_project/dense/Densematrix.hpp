//
// Created by petrov on 12.02.2022.
//

#ifndef SLAE_DENSEMATRIX_HPP
#define SLAE_DENSEMATRIX_HPP
#include <vector>
#include <set>
#include "../utility/Triplet.hpp"
#include "my_project/solvers/ThreeDiadonalSolver.hpp"
#include <iostream>
template<typename T>
class DenseMatrix{
public:
    using elm_t = T;            // Тип данных элементов матрицы
    using idx_t = std::size_t;  // Тип индекса

private:

    std::vector<T> matrix;      //Вектор значений
    idx_t H, W;                 //Размеры матрицы

public:
    /***
     * Конструктор плотной матрицы по ее размерностям, заполнение нулями
     * @param h -- число строк
     * @param w -- число столбцов
     */
    DenseMatrix(const idx_t &h, const idx_t& w): H(h), W(w){
        matrix.resize(H*W);
    };

    /***
     * Конструктор плотной матрицы по ее размерностям и
     * @param h -- число строк
     * @param w -- число столбцов
     * "@param in -- сет триплетов, каждый из которых содержит индексы позиции элемента и его значение
     */
    DenseMatrix(const idx_t &h, const idx_t& w, const std::set<Triplet<T>>& in): H(h), W(w){
        matrix.resize(H * W);
        for(auto k : in){
            matrix[k.i*W + k.j] = k.value;
        }
    };

    /***
     * Оператор получения элемента матрицы по индексам
     * @param i -- Номер строки
     * @param j -- Номер столбца
     * @return Значение элемента в позиции (i, j)
     */
    elm_t& operator()(const idx_t& i, const idx_t& j){
        return matrix[i*W + j];
    };

    /***
     * Оператор получения элемента матрицы по индексам
     * @param i -- Номер строки
     * @param j -- Номер столбца
     * @return Значение элемента в позиции (i, j)
     */
    const elm_t& operator()(const idx_t& i, const idx_t& j) const{
        return matrix[i*W + j];
    };

    /***
     * Метод, позволяющий узнать количество строк матрицы
     * @return количесво строк матрицы
     */
    [[nodiscard]] const idx_t& sizeH() const{
        return H;
    };

    /***
     * Метод, позволяющий узнать количество столбцов матрицы
     * @return количесво столбцов матрицы
     */
    [[nodiscard]] const idx_t& sizeW() const{
        return W;
    };

    /***
     * Метод, меняющий местами две строки матрицы
     * @param first -- первая строчка
     * @param second -- вторая строчка
     */
    void swap(const idx_t& first, const idx_t& second){
        if (H - 1 < first || H - 1 < second) {
            std::stringstream buff;
            buff << "Incorrect indexes of rows given for swap! Rows in matrix: " << H << ". Given first: "
                 << first <<". Second:"<< second << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw Slae::SlaeBaseExceptionCpp(buff.str());
        }
        for(int i = 0; i < W; ++i){
            std::swap(matrix[first * W + i], matrix[second * W + i]);
        }
    };

    /***
     * Метод, удаляющий последнюю строку в матрице
     */
    void deleteLastRow(){
        matrix.erase(matrix.end() - W, matrix.end());
        --H;
    };


};
template<typename T>
std::ostream &operator<<(std::ostream &os, const DenseMatrix<T> &A){
    for(int i = 0; i < A.sizeH(); ++i){
        for(int j = 0; j < A.sizeW(); ++j){
            os << A(i, j) << " ";
        }
        os <<"\n";
    }
    return os;
}
#endif//SLAE_DENSEMATRIX_HPP
