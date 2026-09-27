#include <cstdlib>
#include <ctime>
#include <iostream>

#include "distributions.hpp"
#include "mixture_distribution.hpp"
#include "empiric_distribution.hpp"


void test_3_3_1_main()
{
    const double mu = 1.0;
    const double lambda = 2.0;
    const double nu = 4.0;

    const int sizes[] = {
        200, 1000, 10000, 100000
    };

    std::cout
        << "\n========================================\n"
        << "3.3.1 Основное распределение\n"
        << "========================================\n";

    characteristics theor{};

    theor.mean = mu;
    theor.dispersion =
        dispersion_shift_scale(nu, lambda);
    theor.skewness = 0.0;
    theor.kurtosis =
        coeff_kurtosis(nu);

    print_characteristics(
        "\nТеоретические характеристики",
        theor
    );

    for (int s : sizes) {

        double* sample = new double[s];

        for (int i = 0; i < s; ++i)
            sample[i] =
                modeling_variable_shift_scale(
                    mu, lambda, nu
                );

        characteristics emp =
            empirical_characteristics(
                sample, s
            );

        std::cout
            << "\nВыборка n = "
            << s << "\n";

        print_characteristics(
            "Эмпирические характеристики",
            emp
        );

        print_diff(emp, theor);

        if (s == 1000) {

            empiric ed;

            build_empiric_density(
                &ed,
                sample,
                s
            );

            save_density_points_osn(
                sample,
                s,
                &ed,
                mu,
                lambda,
                nu,
                "density_main.dat"
            );

            generate_theory_curve_main(
                mu,
                lambda,
                nu,
                ed.xl,
                ed.xr,
                500,
                "theory_main.dat"
            );

            plot_density(
                "density_main.dat",
                "theory_main.dat",
                "Основное распределение",
                "density_main.png"
            );

            free_empiric_density(&ed);
        }

        delete[] sample;
    }
}


void test_3_3_1_mixture()
{
    const double mu1 = 0.0;
    const double lambda1 = 1.0;
    const double nu1 = 4.0;

    const double mu2 = 2.0;
    const double lambda2 = 1.0;
    const double nu2 = 2.0;

    const double p = 0.4;

    const int sizes[] = {
        200, 1000, 10000, 100000
    };

    std::cout
        << "\n========================================\n"
        << "3.3.1 Смесь распределений\n"
        << "========================================\n";

    double M1 = mu1;
    double M2 = mu2;

    double D1 =
        dispersion_shift_scale(
            nu1, lambda1
        );

    double D2 =
        dispersion_shift_scale(
            nu2, lambda2
        );

    double gamma11 = 0.0;
    double gamma12 = 0.0;

    double gamma21 =
        coeff_kurtosis(nu1);

    double gamma22 =
        coeff_kurtosis(nu2);

    characteristics theor{};

    theor.mean =
        expectation_mixture(
            M1, M2, p
        );

    theor.dispersion =
        dispersion_mixture(
            M1, D1,
            M2, D2,
            p
        );

    theor.skewness =
        coeff_asymmetry_mixture(
            M1, D1, gamma11,
            M2, D2, gamma12,
            p
        );

    theor.kurtosis =
        coeff_kurtosis_mixture(
            M1, D1, gamma11, gamma21,
            M2, D2, gamma12, gamma22,
            p
        );

    print_characteristics(
        "\nТеоретические характеристики",
        theor
    );

    for (int s : sizes) {

        double* sample = new double[s];

        for (int i = 0; i < s; ++i) {
            sample[i] =
                simulate_mixture(
                    mu1, lambda1, nu1,
                    mu2, lambda2, nu2,
                    p
                );
        }

        characteristics emp =
            empirical_characteristics(
                sample, s
            );

        std::cout
            << "\nВыборка n = "
            << s << "\n";

        print_characteristics(
            "Эмпирические характеристики",
            emp
        );

        print_diff(emp, theor);

        if (s == 1000) {

            empiric ed;

            build_empiric_density(
                &ed,
                sample,
                s
            );

            save_density_points_mix(
                sample,
                s,
                &ed,
                mu1, lambda1, nu1,
                mu2, lambda2, nu2,
                p,
                "density_mixture.dat"
            );

            generate_theory_curve_mixture(
                mu1, lambda1, nu1,
                mu2, lambda2, nu2,
                p,
                ed.xl,
                ed.xr,
                500,
                "theory_mixture.dat"
            );

            plot_density(
                "density_mixture.dat",
                "theory_mixture.dat",
                "Смесь распределений",
                "density_mixture.png"
            );

            free_empiric_density(&ed);
        }

        delete[] sample;
    }
}


void test_3_3_2()
{
    const double mu = 1.0;
    const double lambda = 2.0;
    const double nu = 4.0;

    const int n = 1000;

    std::cout
        << "\n========================================\n"
        << "3.3.2 Эмпирическое распределение\n"
        << "========================================\n";

    characteristics theor{};

    theor.mean = mu;
    theor.dispersion =
        dispersion_shift_scale(
            nu, lambda
        );

    theor.skewness = 0.0;
    theor.kurtosis =
        coeff_kurtosis(nu);

    double* sample =
        new double[n];

    for (int i = 0; i < n; ++i) {
        sample[i] =
            modeling_variable_shift_scale(
                mu, lambda, nu
            );
    }

    empiric ed;

    build_empiric_density(
        &ed,
        sample,
        n
    );

    double* new_sample =
        new double[n];

    modeling_empiric(
        sample,
        n,
        new_sample,
        n
    );

    characteristics original =
        empirical_characteristics(
            sample,
            n
        );

    characteristics empirical =
        empirical_characteristics(
            new_sample,
            n
        );

    print_characteristics(
        "\nТеоретические характеристики",
        theor
    );

    print_characteristics(
        "\nИсходная выборка",
        original
    );

    print_characteristics(
        "\nНовая выборка",
        empirical
    );

    std::cout
        << "\nСравнение новой выборки "
           "с исходной:\n";

    print_diff(
        empirical,
        original
    );

    std::cout
        << "\nСравнение новой выборки "
           "с теоретическими значениями:\n";

    print_diff(
        empirical,
        theor
    );

    free_empiric_density(&ed);

    delete[] sample;
    delete[] new_sample;
}


void tests3()
{
    test_3_3_1_main();
    test_3_3_1_mixture();
    test_3_3_2();
}
