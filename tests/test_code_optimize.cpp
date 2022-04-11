//
// Created by Арсений Плахотнюк on 26.03.2022.
//
#include "gtest/gtest.h"
#include "my_project/sparse/CSR.hpp"
#include <my_project/solvers/CG.hpp>
#include "my_project/utility/Triplet.hpp"
#include<iostream>

#include <array>
#include <string>
#include <iostream>
#include <variant>
#include <iostream>
#include <vector>
#include <functional>


enum COUNTER {
    FIRST = 0,
    SECOND = 1
};


std::array<std::function<void(int)>, 2> functions{
        std::function<void(int)>([](int a){ std::cout << a << std::endl; }),
        std::function<void(int)>([](int a){ std::cout << 2 * a << std::endl; }),
};


void calc(COUNTER counter, int a){
    functions[static_cast<int>(counter)](a);
};


template<int Counter>
struct Calc{
    static void calc(int a);
};

template<>
void Calc<COUNTER::FIRST>::calc(int a){
    std::cout << a << std::endl;
};

template<>
void Calc<COUNTER::SECOND>::calc(int a){
    std::cout << 2 * a << std::endl;
};


TEST(TEST_CODE, NODISCARD){
    calc(COUNTER::SECOND, 5);
    Calc<COUNTER::FIRST>::calc(7);
    Calc<COUNTER::SECOND>::calc(10);
}
