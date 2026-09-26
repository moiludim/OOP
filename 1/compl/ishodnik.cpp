#include <cmath>
#include "spec_func.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>



//функция плотности по варианту
double mainfunc(double x, double nu){
    if(nu<=0){return 0;}
    double result = 1 / (std::pow(2,nu-1) * std::beta(nu/2,nu/2) * std::pow(cosh(x), nu));
    return result;
}

//сдвиг масштабное распределение( как я понимаю это для смеси)
double shift_scale_mainfunc(double x, double mu, double lambda, double nu){
    if(lambda <= 0 or nu <= 0){return 0;}
    double result = (1/lambda) * mainfunc(((x-mu)/lambda), nu);
    return result;
}

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
    



// dispersion
double dispersion(double nu){
    if(nu<=0){
        std::cout<<"nu > 0 ";
        return 0;
    }
    double result = 0.5 * trigamma(nu/2);
    return result;
}

double dispersion_shift_scale(double nu, double lambda){
    if(nu<=0 or lambda <= 0){
        std::cout<<"nu > 0  or lambda > 0 ";
        return 0;
    }
    double result = 0.5 * trigamma(nu/2);
    return result*lambda*lambda;
}

//коэффициент эксцесса
double coeff_kurtosis(double nu){
    double result  = 0.5 * (pentagamma(nu/2)/ (pow(trigamma(nu/2), 2)));
    return result;    
}

/* 

мат ожидание всегда равно 0 и также коэфициент асимметрии = 0

сдвиг масштабное преобразование случайной величины у нас это короче ваще прикол просто рандомим по крутому
случайную величину, почему просто не используем рандом, потому что нам нужны 
именно случайные величины не распределенные нормально

*/

// моделирование с помщью сдвиг масштабного преобразования 
double modeling_variable_shift_scale(double mu, double lambda){

    double r = static_cast<double>(rand()) / RAND_MAX;
    double result = mu + lambda*r;
    return result;
}
//моделирование по варианту случайно величины
double modeling_variable_osnraspr(double nu){
    if(nu<0){std::cout<<"nu > 0";}

    if(nu < 2){
        
        while(true){
            double r1 = static_cast<double>(rand()) / RAND_MAX;
            double x = log(tan(M_PI*r1/2));
            double r2 = static_cast<double>(rand()) / RAND_MAX;
            if(r2 <= std::pow(cosh(x), 1-nu)){
                return x;
            }
            else{
                continue; 
            }
        }

    }

    else{

        while(true){
            double r1 = static_cast<double>(rand()) / RAND_MAX;
            double x = log(0.5*(r1/(1-r1)));
            double r2 = static_cast<double>(rand()) / RAND_MAX;
            if(r2 <= std::pow(cosh(x), 2-nu)){
                return x;
            }
            else{
                continue; 
            }
        }

    }
    
}
//мат ожидание равно смещению для распределния из варика

double expectation_mixture(double M1, double M2, double p){
    if(p < 0 or p > 1){
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

//тесты для основного распределения
void tests1(){

    const double eps = 1e-12;
    std::cout<<"Тесты для основного распределения"<<std::endl;
    std::cout<<"Тест для стандартного распределения при x=0, μ=0, λ=1, ν=4"<<std::endl;
    double test1 = mainfunc(0,4);
    double table_expected_value = 0.750;
    if(abs(test1 - table_expected_value) < eps){
        std::cout<<"Тест пройден!"<<std::endl;
    }
    std::cout<<"Тест для масштабного преобразовани x=0, μ=0, λ=2, ν=4"<<std::endl;
    double test2  = shift_scale_mainfunc(0, 0, 2, 4);
    if(abs(test2 - table_expected_value/2) < eps){
            std::cout<<"Тест пройден!"<<std::endl;
        }
    std::cout<<"Тест для сдвиг-масштабного преобразования распределения при x=μ, μ=5, λ=2, ν=4"<<std::endl;
    double test3 = shift_scale_mainfunc(5, 5, 2, 4);
    if(abs(test3 - table_expected_value/2) < eps){
        std::cout<<"Тест пройден!"<<std::endl;
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



//эмпирическое распределение 

struct empiric {
    double* countn;
    int n; 
    int k;
    double xl, xr;
    double h;
};


// поиск минимума и максимума в массиве
void findMinMax(const double* sample, int n, double* outMin, double* outMax) {
    double mn = sample[0];
    double mx = sample[0];
    for (int j = 1; j < n; ++j) {
        if (sample[j] < mn) mn = sample[j];
        if (sample[j] > mx) mx = sample[j];
    }
    *outMin = mn;
    *outMax = mx;
}

// число промежутков по формуле Стёрджеса: k = 1 + [log2(n)]
int computeSturgesK(int n) {
    return 1 + (int)floor(log2((double)n));
}

// номер промежутка, которому принадлежит x
int binIndex(const empiric* ed, double x) {
    int i = (int)floor((x - ed->xl) / ed->h);
    if (i < 0) i = 0;
    if (i >= ed->k) i = ed->k - 1; // правая граница попадает в последний бин
    return i;
}

// построение структуры EmpiricalDensity по выборке
void buildEmpiricalDensity(empiric* ed, const double* sample, int n, int k_override) {
    if (n <= 0) {
        std::cerr << "buildEmpiricalDensity: empty sample\n";
        exit(1);
    }

    ed->n = n;
    findMinMax(sample, n, &ed->xl, &ed->xr);
    ed->k = (k_override > 0) ? k_override : computeSturgesK(n);
    ed->h = (ed->xr - ed->xl) / ed->k;

    ed->countn = new double[ed->k];
    for (int i = 0; i < ed->k; ++i) ed->countn[i] = 0.0;

    for (int j = 0; j < n; ++j) {
        int i = binIndex(ed, sample[j]);
        ed->countn[i] += 1.0;
    }
}

// освобождение памяти
void freeEmpiricalDensity(empiric* ed) {
    delete[] ed->countn;
    ed->countn = nullptr;
}

// значение эмпирической плотности f*(x)
double empiricalDensityValue(const empiric* ed, double x) {
    if (x < ed->xl or x >= ed->xr) {
        if (x == ed->xr) {
            return ed->countn[ed->k - 1] / (ed->n * ed->h);
        }
        return 0.0;
    }
    int i = binIndex(ed, x);
    return ed->countn[i] / (ed->n * ed->h);
}


void modeling_empiric(double sample[], int n, double new_sample[], int new_n){
    for (int i = 0; i < new_n ; i++){
        int index = rand() % n;
        new_sample[i] = sample[index];
    }
}

int main() {
    std::cout<<"Enter n to make a random sample[n]"<<std::endl;
    int n; 
    std::cin>>n;
    double sample[n];
    double nsample[n+n];
    for(int i = 0; i<n; i++){
        sample[i] = modeling_variable_shift_scale(0,1);
    }
    modeling_empiric(sample, n, nsample, n+n);


    empiric ed;
    buildEmpiricalDensity(&ed, sample, n, 0);

    std::cout << "k = " << ed.k << ", h = " << ed.h << "\n";

    for (double x = ed.xl; x <= ed.xr; x += 0.5) {
        std::cout << "f*(" << x << ") = " << empiricalDensityValue(&ed, x) << "\n";
    }

    freeEmpiricalDensity(&ed);
    return 0;
}

/*

*/





/*
int main(){
    
    srand(static_cast<unsigned int>(time(0)));
    
    for( int i = 0; i<100; i++){
        std::cout<<simulate_mixture(2,4,0.4)<<" ";
    }
    
    tests1();
    tests2();
}
*/