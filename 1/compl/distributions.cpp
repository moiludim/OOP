#include <cmath>
#include <cstdlib>
#include <iostream>

#include "spec_func.hpp"
#include "distributions.hpp"


// Функция плотности
double mainfunc(double x, double nu)
{
    if (nu <= 0) {
        return 0;
    }

    return 1.0 /
        (std::pow(2, nu - 1)
        * std::beta(nu / 2, nu / 2)
        * std::pow(std::cosh(x), nu));
}


// Дисперсия
double dispersion(double nu)
{
    if (nu <= 0) {
        std::cout << "nu > 0\n";
        return 0;
    }

    return 0.5 * trigamma(nu / 2);
}


// Коэффициент эксцесса
double coeff_kurtosis(double nu)
{
    if (nu <= 0) {
        return 0;
    }

    return 0.5 *
        (pentagamma(nu / 2)
        / std::pow(trigamma(nu / 2), 2));
}



double modeling_variable_osnraspr(double nu)
{
    if (nu <= 0) {
        std::cout << "nu > 0\n";
        return 0;
    }

    if (nu < 2) {

        while (true) {

            double r1 =
                static_cast<double>(rand()) / RAND_MAX;

            double x =
                std::log(std::tan(M_PI * r1 / 2));

            double r2 =
                static_cast<double>(rand()) / RAND_MAX;

            if (r2 <= std::pow(std::cosh(x), 1 - nu)) {
                return x;
            }
        }
    }

    else {

        while (true) {

            double r1 =
                static_cast<double>(rand()) / RAND_MAX;

            double x =
                std::log(0.5 * (r1 / (1 - r1)));

            double r2 =
                static_cast<double>(rand()) / RAND_MAX;

            if (r2 <= std::pow(std::cosh(x), 2 - nu)) {
                return x;
            }
        }
    }
}


double shift_scale_mainfunc(
    double x,
    double mu,
    double lambda,
    double nu)
{
    if (lambda <= 0 || nu <= 0) {
        return 0;
    }

    return (1.0 / lambda) *
        mainfunc(
            (x - mu) / lambda,
            nu
        );
}


// Дисперсия сдвиг-масштабного распределения
double dispersion_shift_scale(
    double nu,
    double lambda)
{
    if (nu <= 0 || lambda <= 0) {
        std::cout << "nu > 0 or lambda > 0\n";
        return 0;
    }

    return dispersion(nu) * lambda * lambda;
}


double modeling_variable_shift_scale(
    double mu,
    double lambda)
{
    double r =
        static_cast<double>(rand()) / RAND_MAX;

    return mu + lambda * r;
}