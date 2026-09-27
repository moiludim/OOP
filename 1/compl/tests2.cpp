#include <iostream>
#include "distributions.hpp"
#include "mixture_distribution.hpp"

// 3.2. Минимальный набор тестов для смеси распределений
void tests2()
{
    std::cout << "Введите значение x:";
    double x;
    std::cin >> x;
    std::cout << "\nТесты для смеси распределений\n";

    // 3.2.1 — тривиальный случай: μ1=μ2=x, λ1=λ2=2, ν1=ν2=4, p — любое
    std::cout << "3.2.1 тривиальный случай: mu1=mu2=" << x << ", lambda1=lambda2=2, nu=4, p=0.3\n";
    double p = 0.3;
    double mu = x;
    double lambda = 2, nu = 4;
    double D1 = dispersion_shift_scale(lambda, nu);

    double f = mixture(x, mu, lambda, nu, mu, lambda, nu, p);
    double M = expectation_mixture(mu, mu, p);
    double D = dispersion_mixture(mu, D1, mu, D1, p);

    std::cout << "f = " << f << " (ожидается f(" << x << ") из таблицы для nu=4)\n";
    std::cout << "M = " << M << " (ожидается " << mu << ")\n";
    std::cout << "D = " << D << " (ожидается " << D1 << ")\n\n";

    // 3.2.2 — сдвиговые преобразования: μ1=0, μ2=2, λ1=λ2=1, ν=4, p=0.75
    std::cout << "3.2.2 сдвиговые преобразования: mu1=0, mu2=2, lambda1=lambda2=1, nu=4, p=0.75\n";
    p = 0.75;
    double M1 = 0, M2 = 2;
    D1 = dispersion_shift_scale(1, 4);
    double D2 = dispersion_shift_scale(1, 4);

    f = mixture(0, 0, 1, 4, 2, 1, 4, p);
    M = expectation_mixture(M1, M2, p);
    D = dispersion_mixture(M1, D1, M2, D2, p);

    std::cout << "f = " << f << "\n";
    std::cout << "M = " << M << " (ожидается 1.5)\n";
    std::cout << "D = " << D << " (ожидается 0.75)\n\n";

    // 3.2.3 — масштабные преобразования: μ1=μ2=0, λ1=1, λ2=3, ν=4, p=0.5
    std::cout << "3.2.3 масштабные преобразования: mu1=mu2=0, lambda1=1, lambda2=3, nu=4, p=0.5\n";
    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(1, 4);
    D2 = dispersion_shift_scale(3, 4);

    f = mixture(0, 0, 1, 4, 0, 3, 4, p);
    M = expectation_mixture(M1, M2, p);
    D = dispersion_mixture(M1, D1, M2, D2, p);

    std::cout << "f = " << f << " (ожидается (f1(0)+f1(0)/3)/2)\n";
    std::cout << "M = " << M << " (ожидается 0)\n";
    std::cout << "D = " << D << " (ожидается " << (D1 + D2) / 2 << ")\n\n";

    // 3.2.4 — неравные параметры формы: μ1=μ2=0, λ1=λ2=1, ν1≠ν2, p=0.5
    std::cout << "3.2.4 неравные параметры формы: mu1=mu2=0, lambda1=lambda2=1, nu1=2, nu2=4, p=0.5\n";
    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(1, 2);
    D2 = dispersion_shift_scale(1, 4);

    f = mixture(0, 0, 1, 2, 0, 1, 4, p);
    M = expectation_mixture(M1, M2, p);
    D = dispersion_mixture(M1, D1, M2, D2, p);

    std::cout << "f = " << f << " (ожидается (f1(0)+f2(0))/2)\n";
    std::cout << "M = " << M << " (ожидается 0)\n";
    std::cout << "D = " << D << " (ожидается " << (D1 + D2) / 2 << ")\n\n";
}