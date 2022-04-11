//
// Created by Арсений Плахотнюк on 10.04.2022.
//

#include "gtest/gtest.h"
#include <my_project/SlaeBaseException.hpp>
#include "my_project/solvers/GMRES.hpp"
#include "my_project/dense/Densematrix.hpp"
#include "my_project/utility/Triplet.hpp"

TEST(GMRES, HI) {
    double tolerance = 1e-5;
    std::set<Triplet<double>> data{{0, 0, 7.}, {1, 1, 5.}, {2, 2, 9.}};
    CSR<double> matrix(3, 3, data);
    std::vector<double> col = {0.1, 0.2, 0.3};
    std::vector<double> init = {0., 0., 0.};
    std::vector<double> res = Gmres(matrix, col, init, tolerance);
    ASSERT_NEAR(norm(res - std::vector<double>{0.014, 0.04, 0.033}, NormType::SecondNorm), 0., tolerance);
}