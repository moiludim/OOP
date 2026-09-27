
double mixture(double x, double mu1, double lambda1, double nu1, double mu2, double lambda2, double nu2, double p);
double expectation_mixture(double M1, double M2, double p);
double dispersion_mixture(double M1, double D1, double M2, double D2, double p);
double coeff_asymmetry_mixture(double M1, double D1, double gamma11, double M2, double D2, double gamma12, double p);
double coeff_kurtosis_mixture(double M1, double D1, double gamma11, double gamma21, double M2, double D2, double gamma12, double gamma22, double p);
double simulate_mixture(double mu1, double lambda1, double nu1,double mu2, double lambda2, double nu2,double p);
