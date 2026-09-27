
//  Основное распределение 
double mainfunc(double x, double nu);
double dispersion(double nu);
double coeff_kurtosis(double nu);
double modeling_variable_osnraspr(double nu);
void tests1();

//  Сдвиг-масштабное распределение 
double shift_scale_mainfunc(double x, double mu, double lambda, double nu);
double dispersion_shift_scale(double nu, double lambda);
double modeling_variable_shift_scale(double mu, double lambda, double nu);

//  Смесь распределений 
double mixture(double x, double mu1, double lambda1, double nu1,
                double mu2, double lambda2, double nu2,
                double p);

double expectation_mixture(double M1, double M2, double p);

double dispersion_mixture(
    double M1, double D1,
    double M2, double D2,
    double p);

double coeff_asymmetry_mixture(
    double M1, double D1, double gamma11,
    double M2, double D2, double gamma12,
    double p);

double coeff_kurtosis_mixture(
    double M1, double D1,
    double gamma11, double gamma21,
    double M2, double D2,
    double gamma12, double gamma22,
    double p);

double simulate_mixture(double nu1, double nu2, double p);

void tests2();
