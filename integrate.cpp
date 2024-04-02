#include <iostream>
#include <vector>
#include <functional>

/*** Функция для интегрирования на одном отрезке
* a - начало отезка
* b - конец отрезка
* n - количество точек для квадратуры
* func - функция для интегрирования
***/
[[nodiscard]] double integrateOneSeg(double a, double b, unsigned n, const std::function<double(double)>& func) noexcept 
{
    // P_3 = 1/2 (5 x^3 - 3 x)
    // P_4 = 1/8 (35 x^4 - 30 x^2 +3)
    // P_5 = 1/8 (63 x^5 - 70 x^3 + 15 x)
    std::vector<std::vector<double>> solution = { { sqrt(15.) / (-5.), 0, sqrt(15.) / 5.},
                                                  {(-1.) * sqrt((15. + 2. * sqrt(30.)) / 35.), (-1.) * sqrt((15. - 2. * sqrt(30.)) / 35.), sqrt((15. - 2. * sqrt(30.)) / 35.), sqrt((15. + 2. * sqrt(30.)) / 35.)},
                                                  {(-1.) * sqrt((35. + 2. * sqrt(70.)) / 63.), (-1.) * sqrt((35. - 2. * sqrt(70.)) / 63.), 0., sqrt((35. - 2. * sqrt(70.)) / 63.), sqrt((35. + 2. * sqrt(70.)) / 63.)} };
 
    std::vector<std::vector<double>> w = { {5./9., 8./9., 5./9.},
                                           {(18. - sqrt(30.))/36., (18. + sqrt(30.)) / 36., (18. + sqrt(30.)) / 36., (18. - sqrt(30.)) / 36.},
                                           {(322. + 13.*sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900., 128./225., (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900.} };
    
    double I = 0;
    n = n - 3; 
    for (int i = 0; i < n + 3; i++) {
        I = I + func((b + a) / 2. + (b - a) / 2. * solution[n][i]) * w[n][i];
    }
    I = (b - a) / 2. * I;
    return I;
}

/*** Функция для интегрирования на всем отрезке
* a - начало отезка
* b - конец отрезка
* n - количество точек для квадратуры
* s - количество отрезков
* func - функция для интегрирования
***/
[[nodiscard]] double integrate(double a, double b, unsigned n, unsigned s, const std::function<double(double)>& func) noexcept
{
    double I = 0.;
    for (int i = 0; i < s; i++)
    {
        I = I + integrateOneSeg(a + ((b - a) / s)*i, a + ((b - a)/s) * (i + 1.), n, func);
    }
    return I;
}

double COS(double x) {
    return std::cos(x);
}

int main()
{
    double a = 0, b = 5;
    unsigned n = 3, s;
    double I1 = integrateOneSeg(a, b, n, COS);
    std::cout << I1 << ' ' << sin(b)<< std::endl;

    for (int s = 100; s < 120; s++)
    {
        double I = integrate(a, b, n, s, COS);
        std::cout << "s=" << s << " : " << I << ' ' << sin(b) << std::endl;
    }
}
  