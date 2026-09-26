#include <cmath>
#include "spec_func.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "distributions.hpp"

// смесь 
double mixture(double x, double mu1, double lambda1, double nu1,
                double mu2, double lambda2, double nu2,
                double p)
{
    if(p<0 or p> 1){return 0;}
    double f1 = shift_scale_mainfunc(x, mu1, lambda1, nu1);
    double f2 = shift_scale_mainfunc(x, mu2, lambda2, nu2);

    return (1 - p) * f1 + p * f2;
}

//мат ожидание равно смещению для распределния из варика

double expectation_mixture(double M1, double M2, double p){
    if(p < 0 or p > 1 or (M1 == 0 and M2 == 0)){
        return 0;
    }
    

    return (1-p)*M1 + p*M2;
}
//дисперсия для смеси
double dispersion_mixture(
    double M1, double D1,
    double M2, double D2,
    double p)
{
    if(p < 0 or p > 1){
        return 0;
    }

    double M = expectation_mixture(M1, M2, p);

    return (1-p)*(M1*M1 + D1) + p*(M2*M2 + D2) - M*M;
}

double coeff_asymmetry_mixture(
    double M1, double D1, double gamma11,
    double M2, double D2, double gamma12,
    double p)
{
    if(p < 0 or p > 1){
        return 0;
    }

    double M = expectation_mixture(M1, M2, p);

    double D = dispersion_mixture(
        M1, D1,
        M2, D2,
        p
    );

    double term1 =
        pow(M1 - M, 3)
        + 3*(M1 - M)*D1
        + pow(D1, 1.5)*gamma11;

    double term2 =
        pow(M2 - M, 3)
        + 3*(M2 - M)*D2
        + pow(D2, 1.5)*gamma12;

    return ((1-p)*term1 + p*term2)
           / pow(D, 1.5);
}

double coeff_kurtosis_mixture(
    double M1, double D1,
    double gamma11, double gamma21,
    double M2, double D2,
    double gamma12, double gamma22,
    double p)
{
    if(p < 0 or p > 1){
        return 0;
    }

    double M = expectation_mixture(M1, M2, p);

    double D = dispersion_mixture(
        M1, D1,
        M2, D2,
        p
    );

    double term1 =
        pow(M1 - M, 4)
        + 6*pow(M1 - M, 2)*D1
        + 4*(M1 - M)*pow(D1, 1.5)*gamma11
        + pow(D1, 2)*(gamma21 + 3);

    double term2 =
        pow(M2 - M, 4)
        + 6*pow(M2 - M, 2)*D2
        + 4*(M2 - M)*pow(D2, 1.5)*gamma12
        + pow(D2, 2)*(gamma22 + 3);

    return ((1-p)*term1 + p*term2)
           / pow(D, 2)
           - 3;
}

// моделирование случайной величины для микстуры 

double simulate_mixture(double nu1, double nu2, double p) 
{
    if (p < 0.0 or p > 1.0) { return 0.0; }

    double r = static_cast<double>(rand()) / RAND_MAX;

    if (r < (1.0 - p)) {
        return modeling_variable_osnraspr(nu1);
    } else {
        return modeling_variable_osnraspr(nu2);
    }
}

void tests2(){

    const double eps = 1e-12;
    std::cout<<"Тесты для смеси распределений"<<std::endl;


    // 3.2.1
    double p = 0.3;
    double M1 = 5, M2 = 5;
    double D1 = dispersion_shift_scale(4,2);
    double D2 = dispersion_shift_scale(4,2);

    double test = mixture(5,5,2,4,5,2,4,p);

    if(abs(test - 0.750/2) < eps &&
       abs(expectation_mixture(M1,M2,p) - 5) < eps &&
       abs(dispersion_mixture(M1,D1,M2,D2,p) - D1) < eps){
        std::cout<<"3.2.1 Тест пройден!"<<std::endl;
    }


    // 3.2.2
    p = 0.75;
    M1 = 0; M2 = 2;
    D1 = dispersion_shift_scale(4,1);
    D2 = dispersion_shift_scale(4,1);

    test = mixture(0,0,1,4,2,1,4,p);

    double M = expectation_mixture(M1,M2,p);
    double D = dispersion_mixture(M1,D1,M2,D2,p);

    if(abs(test - (
        0.25*shift_scale_mainfunc(0,0,1,4)
        + 0.75*shift_scale_mainfunc(0,2,1,4)
    )) < eps &&
       abs(M - 1.5) < eps){
        std::cout<<"3.2.2 Тест пройден!"<<std::endl;
    }


    // 3.2.3
    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(4,1);
    D2 = dispersion_shift_scale(4,3);

    test = mixture(0,0,1,4,0,3,4,p);

    M = expectation_mixture(M1,M2,p);
    D = dispersion_mixture(M1,D1,M2,D2,p);

    if(abs(test - (mainfunc(0,4) + mainfunc(0,4)/3)/2) < eps &&
       abs(M) < eps &&
       abs(D - (D1+D2)/2) < eps){
        std::cout<<"3.2.3 Тест пройден!"<<std::endl;
    }


    // 3.2.4
    p = 0.5;
    M1 = 0; M2 = 0;
    D1 = dispersion_shift_scale(2,1);
    D2 = dispersion_shift_scale(4,1);

    test = mixture(0,0,1,2,0,1,4,p);

    M = expectation_mixture(M1,M2,p);
    D = dispersion_mixture(M1,D1,M2,D2,p);

    if(abs(test - (mainfunc(0,2)+mainfunc(0,4))/2) < eps &&
       abs(M) < eps &&
       abs(D - (D1+D2)/2) < eps){
        std::cout<<"3.2.4 Тест пройден!"<<std::endl;
    }
}





