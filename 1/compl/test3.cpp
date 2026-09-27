#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "distributions.hpp"
#include "mixture_distribution.hpp"

// =====================================================================
// 3.3. Тестирование эмпирического распределения и моделирования
//      случайных величин (основное распределение + смесь)
//
// Если функции с такими же именами (build_empiric_density,
// calculate_characteristics и т.п.) уже определены в вашем основном
// файле — удалите дубликаты оттуда, чтобы не было ошибки повторного
// определения при линковке.
// =====================================================================

struct empiric {
    double* countn;
    int n;
    int k;
    double xl, xr;
    double h;
};

struct characteristics {
    double mean;
    double dispersion;
    double skewness;
    double kurtosis;
};

void find_min_max(const double sample[], int n, double* out_min, double* out_max)
{
    double mn = sample[0], mx = sample[0];
    for (int i = 1; i < n; ++i) {
        if (sample[i] < mn) mn = sample[i];
        if (sample[i] > mx) mx = sample[i];
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
    if (i < 0) i = 0;
    if (i >= ed->k) i = ed->k - 1;
    return i;
}

void build_empiric_density(empiric* ed, const double sample[], int n, int k_override = 0)
{
    ed->n = n;
    find_min_max(sample, n, &ed->xl, &ed->xr);
    ed->k = (k_override > 0) ? k_override : sturges_k(n);
    ed->h = (ed->xr - ed->xl) / ed->k;
    ed->countn = new double[ed->k];
    for (int i = 0; i < ed->k; ++i)
        ed->countn[i] = 0.0;
    for (int i = 0; i < n; ++i)
        ed->countn[bin_index(ed, sample[i])]++;
}

void free_empiric_density(empiric* ed)
{
    delete[] ed->countn;
    ed->countn = nullptr;
}

double empiric_density_value(const empiric* ed, double x)
{
    if (x < ed->xl || x > ed->xr) return 0.0;
    int i = (x == ed->xr) ? ed->k - 1 : bin_index(ed, x);
    return ed->countn[i] / (ed->n * ed->h);
}

void modeling_empiric(const double sample[], int n, double new_sample[], int new_n)
{
    for (int i = 0; i < new_n; ++i)
        new_sample[i] = sample[rand() % n];
}

characteristics empirical_characteristics(const double sample[], int n)
{
    characteristics c{};
    for (int i = 0; i < n; ++i)
        c.mean += sample[i];
    c.mean /= n;

    double m2 = 0.0, m3 = 0.0, m4 = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = sample[i] - c.mean;
        m2 += d * d;
        m3 += d * d * d;
        m4 += d * d * d * d;
    }
    c.dispersion = m2 / n;
    if (c.dispersion > 0) {
        c.skewness = m3 / (n * pow(c.dispersion, 1.5));
        c.kurtosis = m4 / (n * c.dispersion * c.dispersion) - 3.0;
    }
    return c;
}

void print_characteristics(const char* label, const characteristics& c)
{
    std::cout << label << ":\n";
    std::cout << "  M  = " << c.mean << "\n";
    std::cout << "  D  = " << c.dispersion << "\n";
    std::cout << "  g1 = " << c.skewness << "\n";
    std::cout << "  g2 = " << c.kurtosis << "\n";
}

void print_diff(const characteristics& emp, const characteristics& theor)
{
    std::cout << "  |M_эмп - M_теор| = " << fabs(emp.mean - theor.mean) << "\n";
    std::cout << "  |D_эмп - D_теор| = " << fabs(emp.dispersion - theor.dispersion) << "\n";
    std::cout << "  |g1_эмп - g1_теор| = " << fabs(emp.skewness - theor.skewness) << "\n";
    std::cout << "  |g2_эмп - g2_теор| = " << fabs(emp.kurtosis - theor.kurtosis) << "\n";
}

// плотности (теоретическая и эмпирическая) в точках выборки —
// основное распределение
void save_density_points_osn(const double sample[], int n, const empiric* ed,
                              double mu, double lambda, double nu,
                              const char* filename)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        std::cerr << "Не удалось создать " << filename << "\n";
        return;
    }
    for (int i = 0; i < n; ++i) {
        double x = sample[i];
        double f_emp = empiric_density_value(ed, x);
        double f_theor = shift_scale_mainfunc(x, mu, lambda, nu);
        fprintf(f, "%lf %lf %lf\n", x, f_emp, f_theor);
    }
    fclose(f);
}

// плотности (теоретическая и эмпирическая) в точках выборки — смесь
void save_density_points_mix(const double sample[], int n, const empiric* ed,
                              double mu1, double lambda1, double nu1,
                              double mu2, double lambda2, double nu2, double p,
                              const char* filename)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        std::cerr << "Не удалось создать " << filename << "\n";
        return;
    }
    for (int i = 0; i < n; ++i) {
        double x = sample[i];
        double f_emp = empiric_density_value(ed, x);
        double f_theor = mixture(x, mu1, lambda1, nu1, mu2, lambda2, nu2, p);
        fprintf(f, "%lf %lf %lf\n", x, f_emp, f_theor);
    }
    fclose(f);
}

// мелкая сетка теоретической плотности (для гладкой кривой на графике,
// в отличие от разреженных точек выборки)
void generate_theory_curve_main(double mu, double lambda, double nu,
                                 double xl, double xr, int steps,
                                 const char* filename)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        std::cerr << "Не удалось создать " << filename << "\n";
        return;
    }
    double margin = 0.1 * (xr - xl);
    xl -= margin;
    xr += margin;
    double step = (xr - xl) / steps;
    for (int i = 0; i <= steps; ++i) {
        double x = xl + i * step;
        fprintf(f, "%lf %lf\n", x, shift_scale_mainfunc(x, mu, lambda, nu));
    }
    fclose(f);
}

void generate_theory_curve_mixture(double mu1, double lambda1, double nu1,
                                    double mu2, double lambda2, double nu2,
                                    double p,
                                    double xl, double xr, int steps,
                                    const char* filename)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        std::cerr << "Не удалось создать " << filename << "\n";
        return;
    }
    double margin = 0.1 * (xr - xl);
    xl -= margin;
    xr += margin;
    double step = (xr - xl) / steps;
    for (int i = 0; i <= steps; ++i) {
        double x = xl + i * step;
        fprintf(f, "%lf %lf\n", x, mixture(x, mu1, lambda1, nu1, mu2, lambda2, nu2, p));
    }
    fclose(f);
}

// строит совмещённый график (теор. кривая + точки эмпирической плотности)
// через gnuplot: сохраняет PNG и сразу открывает окно
void plot_density(const char* points_file, const char* theory_file,
                   const char* title, const char* png_out)
{
    FILE* gp = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent", "w");
    if (gp == nullptr) {
        std::cerr << "Не удалось запустить gnuplot\n";
        return;
    }

    fprintf(gp, "set title '%s'\n", title);
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'f(x)'\n");
    fprintf(gp, "set grid\n");

    fprintf(gp, "set terminal pngcairo size 900,600\n");
    fprintf(gp, "set output '%s'\n", png_out);
    fprintf(gp,
        "plot '%s' using 1:2 with lines lw 2 lc rgb 'blue' title 'Теоретическая плотность', "
        "'%s' using 1:2 with points pt 7 ps 0.6 lc rgb 'red' title 'Эмпирическая плотность (точки выборки)'\n",
        theory_file, points_file);
    fprintf(gp, "unset output\n");

    fprintf(gp, "set terminal wxt\n");
    fprintf(gp,
        "plot '%s' using 1:2 with lines lw 2 lc rgb 'blue' title 'Теоретическая плотность', "
        "'%s' using 1:2 with points pt 7 ps 0.6 lc rgb 'red' title 'Эмпирическая плотность (точки выборки)'\n",
        theory_file, points_file);

    fflush(gp);
    _pclose(gp);
}

// ---------------------------------------------------------------------
// 3.3.1 (а): основное (нестандартное) распределение.
// Несколько выборок разного объёма, сравнение эмпирических и
// теоретических характеристик; для выборки среднего объёма — плотности
// в точках выборки (для построения графика).
// ---------------------------------------------------------------------
void test_3_3_1_main()
{
    const double mu = 1.0, lambda = 2.0, nu = 4.0;
    const int sizes[] = {200, 1000, 10000, 100000};

    std::cout << "\n========================================\n";
    std::cout << "3.3.1 (основное распределение), mu=" << mu
               << " lambda=" << lambda << " nu=" << nu << "\n";
    std::cout << "========================================\n";

    characteristics theor{};
    theor.mean = mu;
    theor.dispersion = dispersion_shift_scale(nu, lambda);
    theor.skewness = 0.0;
    theor.kurtosis = coeff_kurtosis(nu);
    print_characteristics("\nТеоретические характеристики", theor);

    for (int s : sizes) {
        double* sample = new double[s];
        for (int i = 0; i < s; ++i)
            sample[i] = modeling_variable_shift_scale(mu, lambda);

        characteristics emp = empirical_characteristics(sample, s);

        std::cout << "\nВыборка n = " << s << "\n";
        print_characteristics("Эмпирические характеристики", emp);
        print_diff(emp, theor);

        // выборка не слишком большого объёма — считаем плотности в её точках
        if (s == 1000) {
            empiric ed;
            build_empiric_density(&ed, sample, s);
            save_density_points_osn(sample, s, &ed, mu, lambda, nu, "density_main.dat");
            std::cout << "  Плотности (теор./эмп.) в точках выборки -> density_main.dat\n";

            generate_theory_curve_main(mu, lambda, nu, ed.xl, ed.xr, 500, "theory_main.dat");
            plot_density("density_main.dat", "theory_main.dat",
                         "Основное распределение: теория vs эмпирика",
                         "density_main.png");
            std::cout << "  График сохранён в density_main.png\n";

            free_empiric_density(&ed);
        }

        delete[] sample;
    }
}

// ---------------------------------------------------------------------
// 3.3.1 (б): смесь распределений — аналогично основному.
// ---------------------------------------------------------------------
void test_3_3_1_mixture()
{
    const double mu1 = 0.0, lambda1 = 1.0, nu1 = 4.0;
    const double mu2 = 2.0, lambda2 = 1.0, nu2 = 2.0;
    const double p = 0.4;
    const int sizes[] = {200, 1000, 10000, 100000};

    std::cout << "\n========================================\n";
    std::cout << "3.3.1 (смесь распределений), p=" << p << "\n";
    std::cout << "========================================\n";

    double M1 = mu1, M2 = mu2;
    double D1 = dispersion_shift_scale(nu1, lambda1);
    double D2 = dispersion_shift_scale(nu2, lambda2);
    double gamma11 = 0.0, gamma21 = coeff_kurtosis(nu1);
    double gamma12 = 0.0, gamma22 = coeff_kurtosis(nu2);

    characteristics theor{};
    theor.mean = expectation_mixture(M1, M2, p);
    theor.dispersion = dispersion_mixture(M1, D1, M2, D2, p);
    theor.skewness = coeff_asymmetry_mixture(M1, D1, gamma11, M2, D2, gamma12, p);
    theor.kurtosis = coeff_kurtosis_mixture(M1, D1, gamma11, gamma21,
                                             M2, D2, gamma12, gamma22, p);
    print_characteristics("\nТеоретические характеристики", theor);

    for (int s : sizes) {
        double* sample = new double[s];
        for (int i = 0; i < s; ++i)
            sample[i] = simulate_mixture(nu1, nu2, p);

        characteristics emp = empirical_characteristics(sample, s);

        std::cout << "\nВыборка n = " << s << "\n";
        print_characteristics("Эмпирические характеристики", emp);
        print_diff(emp, theor);

        if (s == 1000) {
            empiric ed;
            build_empiric_density(&ed, sample, s);
            save_density_points_mix(sample, s, &ed,
                                     mu1, lambda1, nu1,
                                     mu2, lambda2, nu2, p,
                                     "density_mixture.dat");
            std::cout << "  Плотности (теор./эмп.) в точках выборки -> density_mixture.dat\n";

            generate_theory_curve_mixture(mu1, lambda1, nu1, mu2, lambda2, nu2, p,
                                           ed.xl, ed.xr, 500, "theory_mixture.dat");
            plot_density("density_mixture.dat", "theory_mixture.dat",
                         "Смесь распределений: теория vs эмпирика",
                         "density_mixture.png");
            std::cout << "  График сохранён в density_mixture.png\n";

            free_empiric_density(&ed);
        }

        delete[] sample;
    }
}

// ---------------------------------------------------------------------
// 3.3.2: по эмпирическому распределению одной из выборок моделируется
// новая выборка того же объёма; её характеристики сравниваются с
// характеристиками исходной выборки и с теоретическими.
// ---------------------------------------------------------------------
void test_3_3_2()
{
    const double mu = 1.0, lambda = 2.0, nu = 4.0;
    const int n = 1000;

    std::cout << "\n========================================\n";
    std::cout << "3.3.2 Моделирование по эмпирическому распределению\n";
    std::cout << "========================================\n";

    characteristics theor{};
    theor.mean = mu;
    theor.dispersion = dispersion_shift_scale(nu, lambda);
    theor.skewness = 0.0;
    theor.kurtosis = coeff_kurtosis(nu);

    double* sample = new double[n];
    for (int i = 0; i < n; ++i)
        sample[i] = modeling_variable_shift_scale(mu, lambda);

    empiric ed;
    build_empiric_density(&ed, sample, n);

    double* new_sample = new double[n];
    modeling_empiric(sample, n, new_sample, n);

    characteristics original_emp = empirical_characteristics(sample, n);
    characteristics new_emp = empirical_characteristics(new_sample, n);

    print_characteristics("\nТеоретические характеристики", theor);
    print_characteristics("\nИсходная выборка (эмп. характеристики)", original_emp);
    print_characteristics("\nНовая выборка (по эмпирическому распределению)", new_emp);

    std::cout << "\nСравнение новой выборки с исходной:\n";
    print_diff(new_emp, original_emp);
    std::cout << "\nСравнение новой выборки с теоретическими значениями:\n";
    print_diff(new_emp, theor);

    free_empiric_density(&ed);
    delete[] sample;
    delete[] new_sample;
}

// ---------------------------------------------------------------------
// Точка входа для всех тестов раздела 3.3
// ---------------------------------------------------------------------
void tests3()
{
    test_3_3_1_main();
    test_3_3_1_mixture();
    test_3_3_2();
}

int main(){
    tests3();
}