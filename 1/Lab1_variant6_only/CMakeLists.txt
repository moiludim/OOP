cmake_minimum_required(VERSION 3.15)
project(Lab1_Variant6)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(main_distribution
    main_distribution.cpp
    distributions.cpp
    spec_func.cpp
)

add_executable(mixture_distribution
    mixture_distribution.cpp
    distributions.cpp
    spec_func.cpp
)

add_executable(empiric_distribution
    empiric_distribution.cpp
    mixture_distribution.cpp
    distributions.cpp
    spec_func.cpp
)
