#include <iostream>
#include <cstdlib>

int main()
{
    int choice;

    std::cout << "Что сделать?\n";
    std::cout << "1. Собрать и запустить все тесты (3.1 + 3.2 + 3.3)\n";
    std::cout << "2. Только собрать (без запуска)\n";
    std::cout << "0. Выход\n";
    std::cout << "Ваш выбор: ";

    std::cin >> choice;

    switch (choice) {

        case 1: {
            int build_result = system(
                "g++ main.cpp tests1.cpp tests2.cpp tests3.cpp "
                "mixture_distribution.cpp distributions.cpp empiric_distribution.cpp spec_func.cpp "
                "-o all_tests.exe"
            );

            if (build_result != 0) {
                std::cout << "\nОшибка сборки! Проверьте, что все .cpp файлы лежат рядом.\n";
                break;
            }

            std::cout << "\nСборка успешна, запуск...\n\n";
            system(".\\all_tests.exe");
            break;
        }

        case 2: {
            int build_result = system(
                "g++ main.cpp tests1.cpp tests2.cpp tests3.cpp "
                "mixture_distribution.cpp distributions.cpp spec_func.cpp "
                "-o all_tests.exe"
            );

            if (build_result != 0) {
                std::cout << "\nОшибка сборки!\n";
            } else {
                std::cout << "\nСборка успешна: all_tests.exe создан.\n";
            }
            break;
        }

        case 0:
            std::cout << "Выход.\n";
            break;

        default:
            std::cout << "Неверный выбор!\n";
    }

    return 0;
}