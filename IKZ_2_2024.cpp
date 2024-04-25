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
double* Monotony_condition_2_0(int Jchoice, double* masX, double& alpha, int& step, double eps)
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
    return 0;
}

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
    cout << "Choose method:\n1) Monotony condition\n2) etc\n3) etc\n4) etc\n\nYour choice:";
    cin >> Mchoice;
    separator();
    double alpha; // Шаг для 1-го метода
    int step = 0;
    switch (Mchoice)
    {
    case 1:
        cout << "Choice alpha="; // Переменная равная шагу для поиска.
        cin >> alpha;
        separator();
        Monotony_condition_2_0(Jchoice, masX, alpha, step, eps);
        break;
    case 2:
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