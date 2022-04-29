//
// Created by Арсений Плахотнюк on 10.04.2022.
//

#include "gtest/gtest.h"
#include "my_project/Exceptions/SlaeBaseException.hpp"
#include "my_project/solvers/GMRES.hpp"
#include "my_project/dense/Densematrix.hpp"
#include "my_project/utility/Triplet.hpp"

TEST(GMRES, HI)
{
    double tolerance = 1e-6;
    std::set<Triplet<double>> data{{0, 0, 7.}, {1, 1, 5.}, {2, 2, 9.}};
    CSR<double> matrix(3, 3, data);
    std::vector<double> col = {7., 5., 9.};
    std::vector<double> init = {100., 10., 1000.};
    std::vector<double> res = GMRES<double>(matrix, col, init, tolerance);
    double to_compare = Norm<double, NormType::ThirdNorm>::get_norm(res - std::vector<double>{1., 1., 1.});
    ASSERT_NEAR(to_compare, 0., tolerance);
}
