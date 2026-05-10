cmake_minimum_required(VERSION 3.17)
project(slae)

set(CMAKE_CXX_STANDARD 17)



file(GLOB_RECURSE slea_SOURCES *.h *.cpp *.hpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
if (WIN32)
    set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 -lrt ${CMAKE_CXX_STANDARD_LIBRARIES}")
endif (WIN32)

add_library(slae ${slea_SOURCES} slae/matrix/NonLinearBoundaryValueProblemMatrix5.hpp)
target_link_libraries(slae PUBLIC )

target_include_directories(slae PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
