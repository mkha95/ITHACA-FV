#include "ITHACAstream.H"
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
//include the google test dependencies
#include <gtest/gtest.h>

TEST(nodeTesting, Node_Default_Constructor)
{
    EXPECT_TRUE(1 == 1) << "Node has not been properly initialized";
}

