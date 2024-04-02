#include <iostream> 
#include <functional>
#include <fstream>

double func(double x)
{
    const double pi = 3.14159265358979323846; 
    //return (x - 0.1 * sin(x) - pi / 4);
    //return (tan(x) - 4 * x / pi);
    return (log(cosh(x)));
}

double diff(const std::function<double(double)>& func, double x)
{
    double h = 0.00000001;
    return (func(x + h) - func(x - h)) / (2 * h);
}

[[nodiscard]] double bisectionMethod(double a, double b, const std::function<double(double)>& func, unsigned numberOfIterations) noexcept 
{
    double c = 0.;
    for (int i = 0; i < numberOfIterations; i++)
    {
        c = (a + b) / 2.;
        if (func(c) == 0) return c; 
        else if (func(a) * func(c) < 0) b = c; 
        else if (func(c) * func(b) < 0) a = c;
        else return 9;
    }
    return c;
}

[[nodiscard]] double simpleIterationMethod(double inital, const std::function<double(double)>& func, double tau, unsigned numberOfIterations) noexcept
{
    double x = inital;
    for (int i = 0; i < numberOfIterations; i++) x = x + tau * func(x);
    return x;
}

[[nodiscard]] double newtonMethod(double inital, const std::function<double(double)>& func, unsigned numberOfIterations) noexcept
{
    double x = inital;
    for (int i = 0; i < numberOfIterations; i++) x = x - func(x) / diff(func, x);
    return x;
}

int main()
{
    std::ofstream file;
    file.open("/Users/kiraa/source/repos/non_linear_equation/thirdFunc.txt");
    for (unsigned int i = 0; i <= 100; i++)
    {
        //file << i << ": " << func(bisectionMethod(0, 10, func, i)) << std::endl;
        //file << i << ": " << func(newtonMethod(1, func, i)) << std::endl;
        //file << i << ": " << func(bisectionMethod(0.5, 1.5, func, i)) << std::endl;
        //file << i << ": " << func(simpleIterationMethod(1, func, -0.5, i)) << std::endl;
        //file << i << ": " << func(newtonMethod(1, func, i)) << std::endl;
        file << i << ": " << func(bisectionMethod(0.5, 1.5, func, i)) << std::endl;
        file << i << ": " << func(simpleIterationMethod(1, func, -1, i)) << std::endl;
        file << i << ": " << func(newtonMethod(1, func, i)) << std::endl;
    }
    //for (double j = -3.; j < 0.; j += 0.1)
    //{
    //    file << j << ": " << func(simpleIterationMethod(1, func, j, 10)) << std::endl;
    //}
    return 0;
}
