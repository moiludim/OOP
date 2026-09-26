ЛАБОРАТОРНАЯ РАБОТА №1
Вариант 6

В этой версии сохранена исходная реализация основного распределения
из distributions.cpp без изменений.

Состав:
- main_distribution.cpp — тесты 3.1.1–3.1.3;
- mixture_distribution.cpp/.hpp — смесь и тесты 3.2.1–3.2.4;
- empiric_distribution.cpp — эмпирическое распределение,
  характеристики, моделирование, гистограмма и эксперимент 3.3;
- distributions.cpp/.hpp — исходные функции распределений;
- spec_func.cpp/.hpp — специальные функции;
- build.cpp — выбор программы и сборка.

Сборка через build.cpp:
g++ build.cpp -o build.exe
.\build.exe

Или напрямую:
g++ -std=c++17 main_distribution.cpp distributions.cpp spec_func.cpp -o main_distribution.exe
.\main_distribution.exe

Для смеси:
g++ -std=c++17 mixture_distribution.cpp distributions.cpp spec_func.cpp -o mixture_distribution.exe
.\mixture_distribution.exe

Для эмпирического распределения:
g++ -std=c++17 empiric_distribution.cpp distributions.cpp mixture_distribution.cpp spec_func.cpp -o empiric_distribution.exe
.\empiric_distribution.exe

Для графика нужен установленный gnuplot по пути:
C:\Program Files\gnuplot\bin\gnuplot.exe
