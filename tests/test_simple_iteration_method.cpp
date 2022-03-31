//
// Created by Арсений Плахотнюк on 09.03.2022.
//
#include "gtest/gtest.h"
#include <my_project/solvers/SimpleIteration.hpp>
#include <my_project/utility/Overloads.hpp>
#include </Users/arseniy/Desktop/SLAE4SEM/SLAE4SEM/src/my_project/solvers/steepest_gradient_descent.hpp>
#include </Users/arseniy/Desktop/SLAE4SEM/SLAE4SEM/src/my_project/solvers/CG.hpp>
TEST(SIMPLEITERATION, TEST) {
    double tolerance = 1e-5;
    std::set<Triplet<double>> data{{0, 0, 7.}, {1, 1, 5.}, {2, 2, 9.}};
    CSR<double> my_csr = CSR<double>(3, 3, data);

    std::vector<double> col = {0.1, 0.2, 0.3};
    std::vector<double> init = {0., 0., 0.};

    std::vector<double> res = SimpleIteration(my_csr, col, init, tolerance, 0.15);
    ASSERT_NEAR(norm(res - std::vector<double>{0.014, 0.04, 0.033}, NormType::ThirdNorm), 0., 1e2 * tolerance);
}

TEST(STEEPEST, TEST) {
    double tolerance = 1e-5;
    std::set<Triplet<double>> data{{0, 0, 7.}, {1, 1, 5.}, {2, 2, 9.}};
    CSR<double> my_csr = CSR<double>(3, 3, data);

    std::vector<double> col = {0.1, 0.2, 0.3};
    std::vector<double> init = {0., 0., 0.};
    std::vector<double> answ = {0.014, 0.04, 0.033};

    std::vector<double> res = Steepest_gradient_descent(my_csr, col, init, tolerance);

    ASSERT_NEAR(norm(res - answ, NormType::ThirdNorm), 0., 1e2 * tolerance);
}

TEST(CG, TEST) {
    double tolerance = 1e-5;
    std::set<Triplet<double>> data{{0, 0, 9.}, {1, 1, 9.},
                                   {0, 1, 8.}, {1, 0, 8.}};
    CSR<double> my_csr = CSR<double>(2, 2, data);

    std::vector<double> col = {1., 1.};
    std::vector<double> init = {5., 10.};
//    std::vector<double> answ = {0.0, 0.0};

    std::vector<double> res = CG<double>(my_csr, col, init, tolerance);

//    ASSERT_NEAR(norm(res - answ, NormType::ThirdNorm), 0.,  tolerance);
}