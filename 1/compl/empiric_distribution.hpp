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

void find_min_max(const double sample[], int n,
                  double* out_min, double* out_max);

int sturges_k(int n);

int bin_index(const empiric* ed, double x);

void build_empiric_density(empiric* ed,
                           const double sample[],
                           int n,
                           int k_override = 0);

void free_empiric_density(empiric* ed);

double empiric_density_value(const empiric* ed, double x);

void modeling_empiric(const double sample[],
                      int n,
                      double new_sample[],
                      int new_n);

characteristics empirical_characteristics(const double sample[], int n);

void print_characteristics(const char* label,
                           const characteristics& c);

void print_diff(const characteristics& emp,
                const characteristics& theor);

void save_density_points_osn(
    const double sample[],
    int n,
    const empiric* ed,
    double mu,
    double lambda,
    double nu,
    const char* filename);

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
    const char* filename);

void generate_theory_curve_main(
    double mu,
    double lambda,
    double nu,
    double xl,
    double xr,
    int steps,
    const char* filename);

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
    const char* filename);

void plot_density(
    const char* points_file,
    const char* theory_file,
    const char* title,
    const char* png_out);
