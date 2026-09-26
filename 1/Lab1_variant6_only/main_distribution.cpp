#include <cmath>
#include <iostream>

#include "distributions.hpp"

void tests1()
{
    const double eps = 1e-12;

    std::cout << "Тесты для основного распределения\n";

    double test1 = mainfunc(0, 4);

    if (std::abs(test1 - 0.750) < eps)
        std::cout << "Тест 1 пройден!\n";

    double test2 = shift_scale_mainfunc(0, 0, 2, 4);

    if (std::abs(test2 - 0.750 / 2) < eps)
        std::cout << "Тест 2 пройден!\n";

    double test3 = shift_scale_mainfunc(5, 5, 2, 4);

    if (std::abs(test3 - 0.750 / 2) < eps)
        std::cout << "Тест 3 пройден!\n";
}

int main()
{
    tests1();
    return 0;
}
