
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstdio>
#include "distributions.hpp"
#include "mixture_distribution.hpp"

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

characteristics calculate_characteristics_osn(const double sample[], int n, double nu) {
    characteristics c{};
    c.mean = 0;
    c.skewness = 0;
    c.dispersion = dispersion_shift_scale(nu, 1.0);
    c.kurtosis = coeff_kurtosis(nu);
    return c;
}

characteristics calculate_characteristics(const double sample[], int n) {
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
        c.skewness = m3 / (n * pow(c.dispersion, 1.5));
        c.kurtosis = m4 / (n * c.dispersion * c.dispersion) - 3.0;
    }
    return c;
}

void print_characteristics(const characteristics& c) {
    std::cout << "M  = " << c.mean << "\n";
    std::cout << "D  = " << c.dispersion << "\n";
    std::cout << "g1 = " << c.skewness << "\n";
    std::cout << "g2 = " << c.kurtosis << "\n";
}

void find_min_max(const double sample[], int n, double* out_min, double* out_max) {
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

int sturges_k(int n) {
    return 1 + (int)floor(log2((double)n));
}

int bin_index(const empiric* ed, double x) {
    int i = (int)floor((x - ed->xl) / ed->h);
    if (i < 0)
        i = 0;
    if (i >= ed->k)
        i = ed->k - 1;
    return i;
}

void build_empiric_density(empiric* ed, const double sample[], int n, int k_override) {
    if (n <= 0) {
        std::cerr << "Ошибка: пустая выборка\n";
        exit(1);
    }
    ed->n = n;
    find_min_max(sample, n, &ed->xl, &ed->xr);
    ed->k = (k_override > 0) ? k_override : sturges_k(n);
    ed->h = (ed->xr - ed->xl) / ed->k;
    ed->countn = new double[ed->k];
    for (int i = 0; i < ed->k; ++i)
        ed->countn[i] = 0.0;
    for (int i = 0; i < n; ++i) {
        int bin = bin_index(ed, sample[i]);
        ed->countn[bin]++;
    }
}

void free_empiric_density(empiric* ed) {
    delete[] ed->countn;
    ed->countn = nullptr;
}

double empiric_density_value(const empiric* ed, double x) {
    if (x < ed->xl || x > ed->xr)
        return 0.0;
    if (x == ed->xr)
        return ed->countn[ed->k - 1] / (ed->n * ed->h);
    int i = bin_index(ed, x);
    return ed->countn[i] / (ed->n * ed->h);
}

void modeling_empiric(const double sample[], int n, double new_sample[], int new_n) {
    for (int i = 0; i < new_n; ++i) {
        int index = rand() % n;
        new_sample[i] = sample[index];
    }
}

void print_empiric_density(const empiric* ed) {
    std::cout << "\nЭмпирическая плотность\n";
    std::cout << "---------------------------------------------\n";
    std::cout << "№\tИнтервал\t\tcount\tf*(x)\n";
    std::cout << "---------------------------------------------\n";
    for (int i = 0; i < ed->k; ++i) {
        double left = ed->xl + i * ed->h;
        double right = left + ed->h;
        double density = ed->countn[i] / (ed->n * ed->h);
        std::cout << i + 1 << "\t[" << left << "; " << right;
        if (i == ed->k - 1)
            std::cout << "]";
        else
            std::cout << ")";
        std::cout << "\t" << ed->countn[i] << "\t" << density << "\n";
    }
    std::cout << "---------------------------------------------\n";
}

void experiment_main_distribution() {
    const double mu = 0.0;
    const double lambda = 1.0;
    const double nu = 4.0;
    int sizes[] = {100, 1000, 10000, 100000};
    std::cout << "\n========================================\n";
    std::cout << "3.3.1 Основное распределение\n";
    std::cout << "========================================\n";
    std::cout << "\nТеоретические характеристики:\n";
    double M = mu;
    double D = dispersion_shift_scale(nu, lambda);
    double gamma1 = 0.0;
    double gamma2 = coeff_kurtosis(nu);
    std::cout << "M  = " << M << "\n"
              << "D  = " << D << "\n"
              << "g1 = " << gamma1 << "\n"
              << "g2 = " << gamma2 << "\n";
    for (int s : sizes) {
        double* sample = new double[s];
        for (int i = 0; i < s; ++i) {
            sample[i] = modeling_variable_shift_scale(mu, lambda);
        }
        characteristics c = calculate_characteristics_osn(sample, s, nu);
        std::cout << "\nВыборка n = " << s << "\n";
        print_characteristics(c);
        delete[] sample;
    }
}


void experiment_mixture_distribution() {
    const double mu1 = 0.0, lambda1 = 1.0, nu1 = 4.0;
    const double mu2 = 2.0, lambda2 = 1.0, nu2 = 2.0;
    const double p = 0.4;

    int sizes[] = {100, 1000, 10000, 100000};

    std::cout << "\n========================================\n";
    std::cout << "3.3.1 Распределение смеси\n";
    std::cout << "========================================\n";

    std::cout << "\nТеоретические характеристики:\n";

    double M1 = 0;
    double M2 = 0;
    double D1 = dispersion_shift_scale(nu1, lambda1);
    double D2 = dispersion_shift_scale(nu2, lambda2);
    
    double gamma11 = 0.0, gamma21 = coeff_kurtosis(nu1);
    double gamma12 = 0.0, gamma22 = coeff_kurtosis(nu2);

    double M = expectation_mixture(M1, M2, p);
    double D = dispersion_mixture(M1, D1, M2, D2, p);
    double gamma1 = coeff_asymmetry_mixture(M1, D1, gamma11, M2, D2, gamma12, p);
    double gamma2 = coeff_kurtosis_mixture(M1, D1, gamma11, gamma21, M2, D2, gamma12, gamma22, p);

    std::cout << "M  = " << M << "\n"
              << "D  = " << D << "\n"
              << "g1 = " << gamma1 << "\n"
              << "g2 = " << gamma2 << "\n";

    for (int s : sizes) {
        double* sample = new double[s];

        for (int i = 0; i < s; ++i) {
            sample[i] = simulate_mixture(nu1, nu2, p);
        }

        characteristics c = calculate_characteristics(sample, s);

        std::cout << "\nВыборка n = " << s << "\n";
        print_characteristics(c);

        delete[] sample;
    }
}

void experiment_empiric() {
    const int n = 1000;
    double* sample = new double[n];
    double* new_sample = new double[n];
    for (int i = 0; i < n; ++i) {
        sample[i] = modeling_variable_shift_scale(0, 1);
    }
    modeling_empiric(sample, n, new_sample, n);
    characteristics original = calculate_characteristics(sample, n);
    characteristics empirical = calculate_characteristics(new_sample, n);
    std::cout << "\n========================================\n";
    std::cout << "3.3.2 Эмпирическое распределение\n";
    std::cout << "========================================\n";
    std::cout << "\nИсходная выборка:\n";
    print_characteristics(original);
    std::cout << "\nНовая выборка:\n";
    print_characteristics(empirical);
    delete[] sample;
    delete[] new_sample;
}

void save_density_data(const empiric* ed, const char* filename, double nu) {
    FILE* file = fopen(filename, "w");
    if (file == nullptr) {
        std::cerr << "Ошибка создания файла\n";
        return;
    }
    for (int i = 0; i < ed->k; ++i) {
        double left = ed->xl + i * ed->h;
        double right = left + ed->h;
        double x = (left + right) / 2.0;
        double empirical = ed->countn[i] / (ed->n * ed->h);
        double theoretical = shift_scale_mainfunc(x, 0, 1, nu);
        fprintf(file, "%lf %lf %lf\n", x, empirical, theoretical);
    }
    fclose(file);
}

void plot_empiric_density(const empiric* ed) {
    FILE* data = fopen("histogram.dat", "w");
    if (data == nullptr) {
        std::cerr << "Ошибка создания histogram.dat\n";
        return;
    }
    for (int i = 0; i < ed->k; i++) {
        double left = ed->xl + i * ed->h;
        double right = left + ed->h;
        double center = (left + right) / 2.0;
        double density = ed->countn[i] / (ed->n * ed->h);
        fprintf(data, "%lf %lf\n", center, density);
    }
    fclose(data);
    FILE* gp = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent", "w");
    if (gp == nullptr) {
        std::cerr << "Не удалось запустить gnuplot\n";
        return;
    }
    fprintf(gp, "set title 'Эмпирическая плотность'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'f*(x)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set boxwidth %lf\n", ed->h);
    fprintf(gp, "set style fill solid 0.5\n");
    fprintf(gp, "plot 'histogram.dat' using 1:2 with boxes title 'Эмпирическая плотность'\n");
    fflush(gp);
    _pclose(gp);
}

int main() {
    srand(static_cast<unsigned int>(time(0)));
    
    experiment_main_distribution();
    experiment_empiric();

    const int n = 1000;
    const double nu = 4.0;
    double* sample = new double[n];
    
    for (int i = 0; i < n; ++i) {
        sample[i] = modeling_variable_shift_scale(0, 1);
    }
    
    empiric ed;
    build_empiric_density(&ed, sample, n, 0);
    
    std::cout << "\nk = " << ed.k << "\n";
    std::cout << "h = " << ed.h << "\n";
    std::cout << "xl = " << ed.xl << "\n";
    std::cout << "xr = " << ed.xr << "\n";
    
    print_empiric_density(&ed);
    save_density_data(&ed, "density.dat", nu);
    plot_empiric_density(&ed);
    
    free_empiric_density(&ed);
    delete[] sample;
    
    return 0;
}