#include <iostream>
#include <cmath>

// 6 Вариант

// в т (0,1/2) Значение 1/4 для J1

using namespace std;

// Функции
double J(int choice, double* X) // Функция
{
    switch (choice)
    {
    case 1:
        return 2 * pow(X[0], 2) - 2 * X[0] * X[1] + 3 * pow(X[1], 2) + X[0] - 3 * X[1];
    case 2:
        return 5 * pow(X[0], 2) + 3 * pow(X[1], 2) + 2 * pow(X[2], 2) + 2 * X[0] * X[1] + X[0] * X[2] + X[1] * X[2] + 5 * X[0] + X[2];
    }
}
// Градиент J
double grad_J(int Jchoice, double* X, int component)
{
    switch (Jchoice)
    {
    case 1:
        switch (component)
        {
        case 0:
            return 4 * X[0] - 2 * X[1] + 1;
            break;
        case 1:
            return -2 * X[0] + 6 * X[1] - 3;
            break;
        case 2:
            return 0;
            break;
        }
    case 2:
        switch (component)
        {
        case 0:
            return 10 * X[0] + 2 * X[1] + X[2] + 5;
            break;
        case 1:
            return 6 * X[1] + 2 * X[0] + X[2];
            break;
        case 2:
            return 4 * X[2] + X[0] + X[1] + 1;
            break;
        }
        break;
    }
}
// Метод Из условия монотоности
void Monotony_condition_3_0(int Jchoice, double* masX, double& alpha, int& step, double eps)
{
    double* masY = new double[3];
    while (alpha >= eps)
    {
        bool half_alpha = true;
        masY[0] = masX[0] - alpha * grad_J(Jchoice, masX, 0);
        masY[1] = masX[1] - alpha * grad_J(Jchoice, masX, 1);
        masY[2] = masX[2] - alpha * grad_J(Jchoice, masX, 2);
        if (J(Jchoice, masX) > J(Jchoice, masY))
        {
            for (int i = 0; i < 3; i++)
            {
                masX[i] = masY[i];
            }
            half_alpha = false;
        }
        if (half_alpha)
        {
            alpha /= 2;
        }
        step++;
    }
}
// Метод Хука-Дживса
void Hook_Jeeves(int Jchoice, double* masX, double* masLambda, int& step, double eps)
{
    double* first_masX = new double[3]; // Предыдущая базисная точка
    double* second_masX = new double[3]; // Текущая базисная точка
    double* third_masX = new double[3]; // Следующая базисная точка
    double* previous_masX = new double[3];
    double* current_masX = new double[3];
    bool research = true; // Проверка на получение базисной точки из исследующего поиска
    bool find = true; // Проверка на получение точки из поиска по образцу
    double Lambda = 1;
    for (int i = 0; i < 3; i++) // Задаем начальные условия
    {
        previous_masX[i] = masX[i];
        first_masX[i] = masX[i];
        second_masX[i] = masX[i];
    }
    while (Lambda >= eps)
    {
        // Исследующий поиск
        research = false; // Предполагаем, что исследование ничего не даст
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    if (find) // Если поиск успешен, то ищем нормально около этой точки
                    {
                        current_masX[0] = previous_masX[0] + masLambda[0] * pow(-1, i);
                        current_masX[1] = previous_masX[1] + masLambda[1] * pow(-1, j);
                        current_masX[2] = previous_masX[2] + masLambda[2] * pow(-1, k);
                        if (J(Jchoice, current_masX) < J(Jchoice, second_masX))
                        {
                            for (int swap = 0; swap < 3; swap++)
                            {
                                third_masX[swap] = current_masX[swap]; // Нашли точку, она возможно базисная
                            }
                            research = true; // Исследование дало точку
                        }
                    }
                    else // Если поиск НЕ успешен, то ищем около прошлой базисной точки
                    {
                        current_masX[0] = first_masX[0] + masLambda[0] * pow(-1, i);
                        current_masX[1] = first_masX[1] + masLambda[1] * pow(-1, j);
                        current_masX[2] = first_masX[2] + masLambda[2] * pow(-1, k);
                        if (J(Jchoice, current_masX) < J(Jchoice, second_masX))
                        {
                            for (int swap = 0; swap < 3; swap++)
                            {
                                third_masX[swap] = current_masX[swap]; // Нашли точку, она возможно базисная
                            }
                            research = true; // Исследование дало точку
                        }
                    }
                }
            }
        }
        if (research) // Если нашли базисную точку, то запоминаем ее, а 3-ю забываем.
        {
            for (int i = 0; i < 3; i++)
            {
                first_masX[i] = second_masX[i];
                second_masX[i] = third_masX[i];
            }
        }
        // Поиск по образцу
        find = false; // Предполагаем, что поиск не успешен
        for (int i = 0; i < 3; i++) // Ищем точку по образцу
        {
            current_masX[i] = first_masX[i] + 2 * (second_masX[i] - first_masX[i]);
        }
        if (J(Jchoice, current_masX) < J(Jchoice, second_masX)) // Она меньше текущей базисной?
        {
            for (int i = 0; i < 3; i++)
            {
                previous_masX[i] = current_masX[i]; // Запоминаем эту точку
            }
            find = true; // Поиск дал точку
        }
        Lambda = 0; // Обнуление погрешности
        for (int i = 0; i < 3; i++)
        {
            Lambda += masLambda[i]; // Поиск погрешности
        }
        if (find == false && research == false) // Проверка на уменьшение
        {
            for (int i = 0; i < 3; i++) // Меняем сходимость и новые начальные условия
            {
                masLambda[i] /= 2;
                current_masX[i] = current_masX[i];
                previous_masX[i] = current_masX[i];
                first_masX[i] = current_masX[i];
                second_masX[i] = current_masX[i];
            }
            research = true;
            find = true;
        }
        //cout << fixed;
        //cout.precision(2);
        //cout << research << " " << find << "   ";
        //cout << second_masX[0] << " " << second_masX[1] << " " << second_masX[2] << "    " << masLambda[0] << " " << masLambda[1] << " " << masLambda[2] << endl;
    } 
    for (int i = 0; i < 3; i++) // В массив засовываем ответ
    {
        masX[i] = second_masX[i];
    }
}

void Steepest_descent()
{
}
// Отрез кодов друг от друга
void separator()
{
     cout << "\n==================================================================================\n\n";
}

int main()
{
    int Jchoice, Mchoice;
    double eps = pow(10, -4); 
    double* masX = new double[3];
    cout << "You can choose Eps, but i have eps = 10^-4\n\n";
    separator();
    cout << "Choose J:\n1) 2*x1^2 - 2*x1*x2 + 3*x2^2 + x1 - 3*x2\n2) 5*x1^2 + 3*x2^2 + 2*x3^2 + 2*x1*x2 + x1*x3 + x2*x3 + 5*x1 + x3\n\nYour choice: ";
    cin >> Jchoice;
    separator();
    cout << "Choice x1=";
    cin >> masX[0];
    cout << "Choice x2=";
    cin >> masX[1];
    cout << "Choice x3=";
    cin >> masX[2];
    separator();
    cout << "Choose method:\n1) Monotony condition\n2) etc\n3) Hook-Jeeves\n4) etc\n\nYour choice:";
    cin >> Mchoice;
    separator();
    double alpha; // Шаг для 1-го метода
    int step = 0;
    double* masLambda = new double[3];
    if (Mchoice == 3)
    {
        cout << "Choice Lambda 1=";
        cin >> masLambda[0];
        cout << "Choice Lambda 2=";
        cin >> masLambda[1];
        cout << "Choice Lambda 3=";
        cin >> masLambda[2];
        separator();
    }
    switch (Mchoice)
    {
    case 1:
        cout << "Choice alpha="; // Переменная равная шагу для поиска.
        cin >> alpha;
        separator();
        Monotony_condition_3_0(Jchoice, masX, alpha, step, eps);
        break;
    case 2:
        break;
    case 3:
        Hook_Jeeves(Jchoice, masX, masLambda, step, eps);
        break;
    default:
        cout << "Inccorect choice (method)\n";
        break;
    }
    if (Jchoice != 1)
    {
        cout << "Minimum point x=(" << masX[0] << " ; " << masX[1] << " ; " << masX[2] << ")  J(x)=" << J(Jchoice, masX) << " reached with " << step << " steps\n";
    }
    else
    {
        cout << "Minimum point x=(" << masX[0] << " ; " << masX[1] << " ; Any)" << "  J(x) = " << J(Jchoice, masX) << " reached with " << step << " steps\n";
    }
    return 1;
}