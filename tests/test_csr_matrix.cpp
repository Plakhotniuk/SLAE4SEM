//
// Created by Арсений Плахотнюк on 14.02.2022.
//
#include "gtest/gtest.h"
#include "my_project/sparse/CSR.hpp"
#include "my_project/utility/Triplet.hpp"
#include<iostream>
TEST(TEST_CSR_MATRIX, NODISCARD){
    std::vector<double> v = {5., 6., 7., 9.};
    std::vector<std::size_t> c = {2, 1, 2, 3};
    std::vector<std::size_t> r = {0, 1, 1, 3, 4};
    CSR<double> matrix(4, 5, v, c, r);
    std::cout<< matrix;
    std::set<Triplet<double>> data{{0, 0, 1}, {2, 1, 5}};
    CSR<double> newmatrix(3, 3, data);
    std::cout<< newmatrix;
}

