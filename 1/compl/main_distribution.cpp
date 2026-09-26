#include <cmath>
#include <iostream>

#include "distributions.hpp"


void tests1()
{
    const double eps = 1e-12;
    std::cout << "Введите значение x:";
    double x;
    std::cin>>x;
    std::cout<<"\n Все значения в таблице для x=0";
    std::cout << "Тесты для основного распределения\n";
    std::cout << "3.1.1 тест для стандартного распределения: μ=0, λ=1, ν = 4\n";

    // Стандартное распределение
    double test1 = mainfunc(x, 4);

    std::cout<<"f = "<<test1<<std::endl;

    // Масштабирование
    std::cout << "3.1.2 тест для стандартного распределения: μ=0, λ=2, ν = 4\n";
    double test2 = shift_scale_mainfunc(x, 0, 2, 4);

    std::cout<<"f = "<<test2<<std::endl;
    std::cout<<"Ожидаемое значение это значени из таблицы деленное пополам"<<std::endl;



    // Сдвиг + масштабирование
    std::cout << "3.1.3 тест для стандартного распределения: μ=x, λ=2, ν = 4\n";
    double test3 = shift_scale_mainfunc(x, x, 2, 4);

    std::cout<<"f = "<<test3<<std::endl;
    std::cout<<"Ожидаемое значение это значени из таблицы деленное пополам"<<std::endl;
}


int main()
{
    tests1();
    return 0;
}