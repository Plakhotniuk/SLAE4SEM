//
// Created by Арсений Плахотнюк on 09.02.2022.
//

#include "gtest/gtest.h"
#include "my_project/Exceptions/SlaeBaseException.hpp"
void throwException() { throw Slae::SlaeBaseExceptionCpp("Hi"); }

TEST(EXCEPTION, EXCEPTION_HI) {
    bool isCought = false;
    try
    {
        throwException();
    } catch (const Slae::SlaeBaseExceptionCpp &err)
    {
        isCought = true;
    }
    ASSERT_TRUE(isCought);
}
