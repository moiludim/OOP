#include <cmath>
#include <cstdlib>
#include <iostream>

#include "distributions.hpp"

double mixture(double x, double mu1, double lambda1, double nu1,
               double mu2, double lambda2, double nu2, double p)
{
    if (p < 0 || p > 1) return 0;

    double f1 = shift_scale_mainfunc(x, mu1, lambda1, nu1);
    double f2 = shift_scale_mainfunc(x, mu2, lambda2, nu2);

    return (1 - p) * f1 + p * f2;
}

double expectation_mixture(double M1, double M2, double p)
{
    if (p < 0 || p > 1) return 0;

    return (1 - p) * M1 + p * M2;
}

double dispersion_mixture(
    double M1, double D1,
    double M2, double D2,
    double p)
{
    if (p < 0 || p > 1) return 0;

    double M = expectation_mixture(M1, M2, p);

    return (1 - p) * (M1*M1 + D1)
         + p * (M2*M2 + D2)
         - M*M;
}

double coeff_asymmetry_mixture(
    double M1, double D1, double gamma11,
    double M2, double D2, double gamma12,
    double p)
{
    if (p < 0 || p > 1) return 0;

    double M = expectation_mixture(M1, M2, p);
    double D = dispersion_mixture(M1, D1, M2, D2, p);

    if (D <= 0) return 0;

    double term1 =
        std::pow(M1 - M, 3)
        + 3 * (M1 - M) * D1
        + std::pow(D1, 1.5) * gamma11;

    double term2 =
        std::pow(M2 - M, 3)
        + 3 * (M2 - M) * D2
        + std::pow(D2, 1.5) * gamma12;

    return ((1-p)*term1 + p*term2) / std::pow(D, 1.5);
}

double coeff_kurtosis_mixture(
    double M1, double D1,
    double gamma11, double gamma21,
    double M2, double D2,
    double gamma12, double gamma22,
    double p)
{
    if (p < 0 || p > 1) return 0;

    double M = expectation_mixture(M1, M2, p);
    double D = dispersion_mixture(M1, D1, M2, D2, p);

    if (D <= 0) return 0;

    double term1 =
        std::pow(M1 - M, 4)
        + 6 * std::pow(M1 - M, 2) * D1
        + 4 * (M1 - M) * std::pow(D1, 1.5) * gamma11
        + std::pow(D1, 2) * (gamma21 + 3);

    double term2 =
        std::pow(M2 - M, 4)
        + 6 * std::pow(M2 - M, 2) * D2
        + 4 * (M2 - M) * std::pow(D2, 1.5) * gamma12
        + std::pow(D2, 2) * (gamma22 + 3);

    return ((1-p)*term1 + p*term2) / std::pow(D, 2) - 3;
}

double simulate_mixture(double nu1, double nu2, double p)
{
    if (p < 0 || p > 1) return 0;

    double r = static_cast<double>(rand()) / RAND_MAX;

    if (r < (1 - p))
        return modeling_variable_osnraspr(nu1);

    return modeling_variable_osnraspr(nu2);
}

void tests2()
{
    const double eps = 1e-12;

    std::cout << "Тесты для смеси распределений\n";

    double p = 0.3;
    double M1 = 5, M2 = 5;
    double D1 = dispersion_shift_scale(4, 2);
    double D2 = dispersion_shift_scale(4, 2);

    bool ok =
        std::abs(mixture(5,5,2,4,5,2,4,p) - 0.750/2) < eps &&
        std::abs(expectation_mixture(M1,M2,p) - 5) < eps &&
        std::abs(dispersion_mixture(M1,D1,M2,D2,p) - D1) < eps;

    std::cout << (ok ? "3.2.1 Тест пройден!\n" : "3.2.1 Тест НЕ пройден!\n");

    p = 0.75;
    M1 = 0; M2 = 2;
    D1 = dispersion_shift_scale(4,1);
    D2 = dispersion_shift_scale(4,1);

    double M = expectation_mixture(M1,M2,p);

    ok =
        std::abs(mixture(0,0,1,4,2,1,4,p) -
        (0.25*shift_scale_mainfunc(0,0,1,4) +
         0.75*shift_scale_mainfunc(0,2,1,4))) < eps &&
        std::abs(M - 1.5) < eps;

    std::cout << (ok ? "3.2.2 Тест пройден!\n" : "3.2.2 Тест НЕ пройден!\n");

    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(4,1);
    D2 = dispersion_shift_scale(4,3);

    M = expectation_mixture(M1,M2,p);
    double D = dispersion_mixture(M1,D1,M2,D2,p);

    ok =
        std::abs(mixture(0,0,1,4,0,3,4,p) -
        (mainfunc(0,4) + mainfunc(0,4)/3)/2) < eps &&
        std::abs(M) < eps &&
        std::abs(D - (D1+D2)/2) < eps;

    std::cout << (ok ? "3.2.3 Тест пройден!\n" : "3.2.3 Тест НЕ пройден!\n");

    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(2,1);
    D2 = dispersion_shift_scale(4,1);

    M = expectation_mixture(M1,M2,p);
    D = dispersion_mixture(M1,D1,M2,D2,p);

    ok =
        std::abs(mixture(0,0,1,2,0,1,4,p) -
        (mainfunc(0,2) + mainfunc(0,4))/2) < eps &&
        std::abs(M) < eps &&
        std::abs(D - (D1+D2)/2) < eps;

    std::cout << (ok ? "3.2.4 Тест пройден!\n" : "3.2.4 Тест НЕ пройден!\n");
}
