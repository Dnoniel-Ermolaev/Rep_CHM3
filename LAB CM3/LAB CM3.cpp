// Copyright Group 3  2022

#include <iostream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <typeinfo>
#include <cmath>

double equation(double x, double y) { return ((x * x * x) + (y * y * y) + 3.0); }
double deltaU(double x, double y) { return -6.0 * x - 6.0 * y; }

double GU1(double x, double y,double a) {
    return(y * y * y + a*a*a+ 3.0);
}
double GU2(double x, double y,double b) {
    return(y * y * y + b * b * b + 3.0);
}
double GU3(double x, double y,double c) {
    return(x * x * x +c*c*c + 3.0);
}
double GU4(double x, double y,double d) {
    return (x * x * x +d*d*d+ 3.0);
}

std::vector<double> startValueX(double m, double a, double b) {
    double h = (b - a) / m;
    std::vector<double> xi;
    for (int i = 0; i < m + 1; i++) {
        xi.push_back(0 + i * h);
    }
    return xi;
}
std::vector<double> startValueY(double n, double c, double d) {
    double k = (d - c) / n;
    std::vector<double> yi;
    for (int i = 0; i < n + 1; i++) {
        yi.push_back(0 + i * k);
    }
    return yi;
}

void print1(std::vector<std::vector<double> > M, double n, double m) {
    std::cout << std::endl
        << "__________________________________________________________";
    std::cout << std::endl << std::endl;
    for (int i = 0; i < n + 1; i++) {
        // std::cout << "x" << i<<"|";
        std::cout << std::setw(6) << "y" << abs(n - i) << "|";
        for (int j = 0; j < m + 1; j++) {
            std::cout << std::setw(9) << M[j][n - i] << "|";
        }
        std::cout << std::endl;
        std::cout << "__________________________________________________________";
        std::cout << std::endl;
    }
    std::cout << std::setw(8) << "|";
    for (int i = 0; i < m + 1; i++) {
        std::cout << std::setw(8) << "x" << abs(i) << "|";
    }
    std::cout << "\n";
}

void print(std::vector<std::vector<double> > M, double n, double m) {
    for (int i = n; i >= 0; i--) {
        std::cout << std::setw(8) << "x" << abs(i - n) << "|";
    }
    std::cout << std::endl
        << "__________________________________________________________";
    std::cout << std::endl << std::endl;
    for (int j = 0; j < m + 1; j++) {
        // std::cout << "x" << i<<"|";
        for (int i = 0; i < n + 1; i++) {
            std::cout << std::setw(9) << M[i][j] << "|";
        }
        std::cout << std::endl;
        std::cout << "__________________________________________________________";
        std::cout << std::endl;
    }
}

std::vector<std::vector<double> > TrueSolution(std::vector<std::vector<double> > M,
    int n, int m, double a, double b,
    double c, double d) {
    double h = (b - a) / n;
    double k = (d - c) / m;
    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            double xi = a + i * h;
            double yi = c + j * k;
                M[i][j] = equation(xi, yi);
        }
    }
    return M;
}
std::vector<std::vector<double> > rightPart(std::vector<std::vector<double> > F,
    int n, int m, double a, double b,
    double c, double d) {
    double h = (b - a) / n;
    double k = (d - c) / m;
    double h2 = ((n / (b - a)) * (n / (b - a)));
    double k2 = ((m / (d - c)) * (m / (d - c)));

    double a2 = -2.0 * (h2 + k2);

    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            double xi, yi, sum = 0.0;
            xi = a + i * h;
            yi = c + j * k;
            if (j == 1)
                sum += (1.0 / (k * k)) * GU3(xi, yi,c);
            else
                if (j == m - 1)
                    sum += (1.0 / (k * k)) * GU4(xi, yi,d);
            if (i == 1)
                sum += (1.0 / (h * h)) * GU1(xi, yi,a);
            else
                if (i == n - 1)
                    sum += (1.0 / (h * h)) * GU2(xi, yi,b);
            F[i][j] = -deltaU(xi, yi) - sum;
        }
    }

    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            std::cout << std::setw(6) << "F[" << i << "][" << j << "]" << F[i][j]
                << "\n";
        }
    }
    return F;
}
std::vector<std::vector<double> > firstApproc(
    std::vector<std::vector<double> > V, int n, int m, double a, double b,
    double c, double d) {
    double h = (b - a) / n;
    double k = (d - c) / m;
    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            double xi = a + i * h;
            double yi = c + j * k;
            if (i == 0 || j == 0 || i == n || j == m) {
                    V[i][j] = equation(xi, yi);
            }
            else
                V[i][j] = 0;
        }
    }
    return V;
}
double discrepancyCRACK(std::vector<std::vector<double> > M, int n, int m,
    double a, double b, double c, double d) {
    double ep = 0.0;  
    double h = (b - a) / n;
    double k = (d - c) / m;
    double h2 = ((n / (b - a)) * (n / (b - a)));
    double k2 = ((m / (d - c)) * (m / (d - c)));
    double a2 = -2.0 * (h2 + k2);

    std::vector<std::vector<double> > F(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            F[i].push_back(0);
        }
    }
    F = rightPart(F, n, m, a, b,c,d);

    long int n2 = (n - 1) * (m - 1);
    double** A = new double* [n2];
    for (int i = 0; i < n2; i++) {
        A[i] = new double[n2];

        for (int j = 0; j < n2; j++) {
            A[i][j] = 0;
        }
    }

    for (int i = 0; i < (n - 1) * (m - 1); i++) {
        int J1 = i % (n - 1) + 1;
        int I1 = i / (m - 1) + 1;
        A[i][i] = -2 * (h2 + k2);
        if (I1 > 1) A[i][i - n + 1] = k2;
        if (I1 < n - 1) A[i][i + n - 1] = k2;
        if (J1 > 1) A[i][i - 1] = h2;
        if (J1 < n - 1) A[i][i + 1] = h2;
    }

    double r = 0;
    std::vector<double> V;
    std::vector<double> Fr;

    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            V.push_back(M[i][j]);
            Fr.push_back(F[i][j]);
        }
    }
    for (int j = 0; j < n2; j++) {
        // std::cout << V[j]<<"\t\n";;
    }

    for (int i = 0; i < (n - 1) * (m - 1); i++) {
        double temp = 0;

        for (int j = 0; j < (n - 1) * (m - 1); j++) {
            temp += A[i][j] * V[j];
        }
        //std::cout << "AV[" << i << "]=" << abs(temp) << "\tR[" << i
          //  << "]=" << abs(temp - Fr[i]) << "\n";

        r += abs(temp - Fr[i]) * abs(temp - Fr[i]);
    }

    F.clear();
    return sqrt(r);
}

std::vector<std::vector<double> > Zeidal(std::vector<std::vector<double> > M,
    int n, int m, double a, double b,
    double c, double d, double eps,
    int nMax, double& epsMax, int& iter) {
    double curEps = 0;
    double h = (b - a) / (double)n;
    double k = (d - c) / (double)m;
    double h2 = -(((double)n / (b - a)) * ((double)n / (b - a)));
    double k2 = -(((double)m / (d - c)) * ((double)m / (d - c)));
    double a2 = -2.0 * (h2 + k2);

    std::vector<double> xi = startValueX(n, a, b);
    std::vector<double> yi = startValueY(m, c, d);
    double vOld, vCurr;
    double R = 0;
    iter = 0;
    for (int Nn = 0; Nn < nMax; ++Nn) {
            epsMax = 0.0;
            for (int j = 1; j < m; j++) {
                for (int i = 1; i < n; i++) {
                    vOld = M[i][j];

                    vCurr = -(h2 * (M[i + 1][j] + M[i - 1][j]) + k2 * (M[i][j + 1] + M[i][j - 1]));
                    vCurr = vCurr + deltaU(xi[i], yi[i]);
                    vCurr = vCurr / a2;

                    curEps = std::max(epsMax, abs(vOld - vCurr));
                    epsMax = curEps;
                    //std::cout << "\nEPS for[ " << iter << "]=" << curEps << std::endl;

                    M[i][j] = vCurr;
                    //print1(M, n, m);
                }
            }
            iter++;
            if ((epsMax < eps)) break;
    }
    return M;
}

double EpsSLAU(std::vector<std::vector<double> > M, int n, int m, double a,
    double b, double c, double d) {
    std::vector<std::vector<double> > M1(n + 1);
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            M1[i].push_back(0);
        }
    }
    double zs = 0;
    double h = (b - a) / n;
    double k = (d - c) / m;

    M1 = TrueSolution(M1, n, m, a, b, c, d);

    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            double z = abs(M1[i][j] - M[i][j]);
            if (z > zs) zs = z;
            //std::cout << "Z=" << z << "\n";
        }
    }
    M1.clear();
    return zs;
}

int main() {
    setlocale(LC_ALL, "rus");

    int N_max = 1000;   // Максимальное число итераций
    int iter = 0;       // Счетчик итераций
    double eps=1e-6;         // Параметр требуемой точности
    double epsMax = 0;  // Достигнутая точность
    double n = 4;
    double m = 4;
    // std::vector<std::vector<double> > M;  // Искомый вектор
    double a = 0, b = 1.0, c = 0, d = 0.02;  // гранциы
    double disp = 0;
    int FLAG;
    int FLAG2;
    std::cout << "Если желаете ввести значения самостоятельно, введите цифру 1, иначе введите любой символ" << std::endl;
    std::cin >> FLAG;
    if (FLAG == 1) {
        std::cout << "введите Число итераций N: " << std::endl;
        std::cin >> N_max;
        std::cout << "введите границу сетки A: " << std::endl;
        std::cin >> a;
        std::cout << "введите границу сетки B: " << std::endl;
        std::cin >> b;
        std::cout << "введите границу сетки C: " << std::endl;
        std::cin >> c;
        std::cout << "введите границу сетки D: " << std::endl;
        std::cin >> d;
        std::cout << "введите размер сетки N: " << std::endl;
        std::cin >> n;
        std::cout << "введите размер сетки M: " << std::endl;
        std::cin >> m;
        std::cout << std::endl;
        std::cout << "введите параметр требуемой точности: " << std::endl;
        std::cin >> eps;
        std::cout << std::endl;
    }
    
    std::cout << "Желаете выводить таблицы? Да - 1; Нет- любой другой символ" << std::endl;
    std::cin >> FLAG2;

    std::cout
        << "Точное решение задачи Дирихле для уравнения Пуассона во всех узлах: "
        << std::endl;

    std::vector<std::vector<double> > M(n + 1);
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            M[i].push_back(0);
        }
    }
    M = TrueSolution(M, n, m, a, b, c, d);
    if (FLAG2 == 1) {
        if (n == m)
            print1(M, n, m);
        else
            print(M, n, m);
    }
    std::vector<std::vector<double> > FA(n + 1);
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            FA[i].push_back(0);
        }
    }
    std::cout << std::endl;
    std::cout << "Граничные условия и первое приближение: " << std::endl;
    FA = firstApproc(FA, n, m, a, b, c, d);
    if (FLAG2 == 1) {
        if (n == m)
            print1(FA, n, m);
        else
            print(FA, n, m);
    }
    std::cout << std::endl << std::endl;
    std::cout << "Применение метода Зейдаля: " << std::endl;
    FA = Zeidal(FA, n, m, a, b, c, d, eps, N_max, epsMax, iter);
    std::cout << std::endl;
    std::cout << "Последняя итерация метода при заданной точности: (eps = " << eps
        << " )" << std::endl;
    if (FLAG2 == 1) {
        if (n == m)
            print1(FA, n, m);
        else
            print(FA, n, m);
    }

    std::vector<std::vector<double> > F(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            F[i].push_back(0);
        }
    }
    F = rightPart(F, n, m, a, b, c, d);

    std::cout << "================" << std::endl;
    std::cout << "    Результат   " << std::endl;
    std::cout << "================" << std::endl;
    std::cout << "Исходное уравнение: " << std::endl;
    std::cout << "u = x^3+y^3+2" << std::endl;
    std::cout << "Значения по Х: от " << a << " до " << b << std::endl;
    std::cout << "Значения по Y: от " << c << " до " << d << std::endl;
    std::cout << "Размерность сетки: " << n << " x " << m << std::endl;
    std::cout << "Шаг сетки по оси x: h = " << (b - a) / n << std::endl;
    std::cout << "Шаг сетки по оси y: k = " << (d - c) / m << std::endl
        << std::endl;

    std::cout << "Критерии остановки счета:" << std::endl;
    std::cout << "Максимальное число шагов: " << N_max << std::endl;
    std::cout << "Заданная точность: " << eps << std::endl << std::endl;

    std::cout << "Полученные результаты: " << std::endl;
    std::cout << "Число затраченных шагов: " << iter << std::endl;
    std::cout << "Точность на выходе: " << epsMax << std::endl;
    disp = discrepancyCRACK(FA, n, m, a, b, c, d);
    std::cout << "Невязка на выходе (евкл.): " << disp << std::endl;
    double err = EpsSLAU(FA, n, m, a, b, c, d);
    std::cout << "Погрешность решения СЛАУ: " << err << std::endl;
    FA.clear();
    return 0;
}
