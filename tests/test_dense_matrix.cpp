#include "gtest/gtest.h"
#include "my_project/sparse/CSR.hpp"
#include "my_project/utility/Triplet.hpp"
#include "my_project/dense/Densematrix.hpp"
#include<iostream>
TEST(TEST_DENSE_MATRIX, NODISCARD){
    DenseMatrix<double> zero_matrix(5, 5);
    std::cout<< zero_matrix;

    std::cout<<"----"<<std::endl;
    std::set<Triplet<double>> data{{2, 0, 7}, {1, 1, 5}, {0, 3, 9}};
    DenseMatrix<double> matrix(4, 4, data);
    std::cout<< matrix;

    std::cout<<"----"<<std::endl;
    matrix.deleteLastRow();
    std::cout<< matrix;

    std::cout<<"----"<<std::endl;
    matrix.swap(0, 2);
    std::cout<< matrix;
}

