#include <iostream>
#include <cstdlib>

int main()
{
    int choice;

    std::cout << "Выберите программу:\n";
    std::cout << "1. main_distribution.cpp\n";
    std::cout << "2. mixture_test.cpp\n";
    std::cout << "3. empiric_distribution.cpp\n";
    std::cout << "Ваш выбор: ";

    std::cin >> choice;

    switch (choice) {

        case 1:
            system(
                "g++ main_distribution.cpp distributions.cpp spec_func.cpp "
                "-o main_distribution.exe"
            );

            system(".\\main_distribution.exe");
            break;


        case 2:
            system(
                "g++ mixture_test.cpp mixture_distribution.cpp distributions.cpp spec_func.cpp "
                "-o mixture_distribution.exe"
            );

            system(".\\mixture_distribution.exe");
            break;


        case 3:
            system(
                "g++ empiric_distribution.cpp distributions.cpp mixture_distribution.cpp spec_func.cpp -o empiric_distribution.exe "
            );

            system(".\\empiric_distribution.exe");
            break;


        default:
            std::cout << "Неверный выбор!\n";
    }

    return 0;
}