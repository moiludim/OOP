#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "distributions.hpp"
#include "mixture_distribution.hpp"
#include "empiric_distribution.hpp"


void find_min_max(const double sample[], int n,
                  double* out_min, double* out_max)
{
    double mn = sample[0];
    double mx = sample[0];

    for (int i = 1; i < n; ++i) {
        if (sample[i] < mn)
            mn = sample[i];

        if (sample[i] > mx)
            mx = sample[i];
    }

    *out_min = mn;
    *out_max = mx;
}


int sturges_k(int n)
{
    return 1 + (int)floor(log2((double)n));
}


int bin_index(const empiric* ed, double x)
{
    int i = (int)floor((x - ed->xl) / ed->h);

    if (i < 0)
        i = 0;

    if (i >= ed->k)
        i = ed->k - 1;

    return i;
}


void build_empiric_density(empiric* ed,
                           const double sample[],
                           int n,
                           int k_override)
{
    if (n <= 0) {
        std::cerr << "Ошибка: пустая выборка\n";
        exit(1);
    }

    ed->n = n;

    find_min_max(sample, n, &ed->xl, &ed->xr);

    ed->k = (k_override > 0)
                ? k_override
                : sturges_k(n);

    ed->h = (ed->xr - ed->xl) / ed->k;

    ed->countn = new double[ed->k];

    for (int i = 0; i < ed->k; ++i)
        ed->countn[i] = 0.0;

    for (int i = 0; i < n; ++i) {
        int bin = bin_index(ed, sample[i]);
        ed->countn[bin]++;
    }
}


void free_empiric_density(empiric* ed)
{
    delete[] ed->countn;
    ed->countn = nullptr;
}


double empiric_density_value(const empiric* ed, double x)
{
    if (x < ed->xl || x > ed->xr)
        return 0.0;

    int i;

    if (x == ed->xr)
        i = ed->k - 1;
    else
        i = bin_index(ed, x);

    return ed->countn[i] / (ed->n * ed->h);
}


void modeling_empiric(const double sample[],
                      int n,
                      double new_sample[],
                      int new_n)
{
    for (int i = 0; i < new_n; ++i)
        new_sample[i] = sample[rand() % n];
}


characteristics empirical_characteristics(
    const double sample[],
    int n)
{
    characteristics c{};

    for (int i = 0; i < n; ++i)
        c.mean += sample[i];

    c.mean /= n;

    double m2 = 0.0;
    double m3 = 0.0;
    double m4 = 0.0;

    for (int i = 0; i < n; ++i) {
        double d = sample[i] - c.mean;

        m2 += d * d;
        m3 += d * d * d;
        m4 += d * d * d * d;
    }

    c.dispersion = m2 / n;

    if (c.dispersion > 0) {
        c.skewness =
            m3 / (n * pow(c.dispersion, 1.5));

        c.kurtosis =
            m4 / (n * c.dispersion * c.dispersion) - 3.0;
    }

    return c;
}


void print_characteristics(
    const char* label,
    const characteristics& c)
{
    std::cout << label << ":\n";
    std::cout << "  M  = " << c.mean << "\n";
    std::cout << "  D  = " << c.dispersion << "\n";
    std::cout << "  g1 = " << c.skewness << "\n";
    std::cout << "  g2 = " << c.kurtosis << "\n";
}


void print_diff(const characteristics& emp,
                const characteristics& theor)
{
    std::cout
        << "  |M_эмп - M_теор| = "
        << fabs(emp.mean - theor.mean) << "\n";

    std::cout
        << "  |D_эмп - D_теор| = "
        << fabs(emp.dispersion - theor.dispersion) << "\n";

    std::cout
        << "  |g1_эмп - g1_теор| = "
        << fabs(emp.skewness - theor.skewness) << "\n";

    std::cout
        << "  |g2_эмп - g2_теор| = "
        << fabs(emp.kurtosis - theor.kurtosis) << "\n";
}


void save_density_points_osn(
    const double sample[],
    int n,
    const empiric* ed,
    double mu,
    double lambda,
    double nu,
    const char* filename)
{
    FILE* f = fopen(filename, "w");

    if (!f) {
        std::cerr << "Не удалось создать "
                  << filename << "\n";
        return;
    }

    for (int i = 0; i < n; ++i) {
        double x = sample[i];

        double f_emp =
            empiric_density_value(ed, x);

        double f_theor =
            shift_scale_mainfunc(
                x, mu, lambda, nu);

        fprintf(f, "%lf %lf %lf\n",
                x, f_emp, f_theor);
    }

    fclose(f);
}


void save_density_points_mix(
    const double sample[],
    int n,
    const empiric* ed,
    double mu1,
    double lambda1,
    double nu1,
    double mu2,
    double lambda2,
    double nu2,
    double p,
    const char* filename)
{
    FILE* f = fopen(filename, "w");

    if (!f) {
        std::cerr << "Не удалось создать "
                  << filename << "\n";
        return;
    }

    for (int i = 0; i < n; ++i) {
        double x = sample[i];

        double f_emp =
            empiric_density_value(ed, x);

        double f_theor =
            mixture(
                x,
                mu1, lambda1, nu1,
                mu2, lambda2, nu2,
                p);

        fprintf(f, "%lf %lf %lf\n",
                x, f_emp, f_theor);
    }

    fclose(f);
}


void generate_theory_curve_main(
    double mu,
    double lambda,
    double nu,
    double xl,
    double xr,
    int steps,
    const char* filename)
{
    FILE* f = fopen(filename, "w");

    if (!f) {
        std::cerr << "Не удалось создать "
                  << filename << "\n";
        return;
    }

    double margin = 0.1 * (xr - xl);

    xl -= margin;
    xr += margin;

    double step = (xr - xl) / steps;

    for (int i = 0; i <= steps; ++i) {
        double x = xl + i * step;

        fprintf(
            f,
            "%lf %lf\n",
            x,
            shift_scale_mainfunc(
                x, mu, lambda, nu)
        );
    }

    fclose(f);
}


void generate_theory_curve_mixture(
    double mu1,
    double lambda1,
    double nu1,
    double mu2,
    double lambda2,
    double nu2,
    double p,
    double xl,
    double xr,
    int steps,
    const char* filename)
{
    FILE* f = fopen(filename, "w");

    if (!f) {
        std::cerr << "Не удалось создать "
                  << filename << "\n";
        return;
    }

    double margin = 0.1 * (xr - xl);

    xl -= margin;
    xr += margin;

    double step = (xr - xl) / steps;

    for (int i = 0; i <= steps; ++i) {
        double x = xl + i * step;

        fprintf(
            f,
            "%lf %lf\n",
            x,
            mixture(
                x,
                mu1, lambda1, nu1,
                mu2, lambda2, nu2,
                p)
        );
    }

    fclose(f);
}


void plot_density(
    const char* points_file,
    const char* theory_file,
    const char* title,
    const char* png_out)
{
    FILE* gp =
        _popen(
            "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent",
            "w"
        );

    if (gp == nullptr) {
        std::cerr << "Не удалось запустить gnuplot\n";
        return;
    }

    fprintf(gp, "set title '%s'\n", title);
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'f(x)'\n");
    fprintf(gp, "set grid\n");

    fprintf(gp,
        "set terminal pngcairo size 900,600\n");

    fprintf(gp,
        "set output '%s'\n",
        png_out);

    fprintf(
        gp,
        "plot '%s' using 1:2 with lines lw 2 "
        "title 'Теоретическая плотность', "
        "'%s' using 1:2 with points pt 7 "
        "title 'Эмпирическая плотность'\n",
        theory_file,
        points_file
    );

    fprintf(gp, "unset output\n");

    fprintf(gp, "set terminal wxt\n");

    fprintf(
        gp,
        "plot '%s' using 1:2 with lines lw 2 "
        "title 'Теоретическая плотность', "
        "'%s' using 1:2 with points pt 7 "
        "title 'Эмпирическая плотность'\n",
        theory_file,
        points_file
    );

    fflush(gp);
    _pclose(gp);
}