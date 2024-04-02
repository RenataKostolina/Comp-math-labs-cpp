#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

class NewtonianInterpolator{
private:
    double** matrix;
    int size;

public:
     explicit NewtonianInterpolator(const std::vector<double>& x, const std::vector<double>& y) : size(int(x.size())) {
         matrix = new double* [size];
         for (int i = 0; i < size; i++) {
             matrix[i] = new double[i + 2] ;
             matrix[i][0] = x[i];
             matrix[i][1] = y[i];
         }
         for (int i = 2; i < size+2; i++) {
             for (int j = i - 1; j < size; j++) {
                 double p = (matrix[j][i-1] - matrix[j-1][i - 1]);
                 double q = (matrix[j][0] - matrix[j - i + 1][0]);
                 matrix[j][i] = p / q;
             }
         }
     }
     [[nodiscard]] double interpolate(double x) const {
         double func = matrix[size-1][size];
         for (int i = 1; i < size; i++) {
             double p = matrix[size - i - 1][size - i];
             func = func * (x - matrix[size - i - 1][0]) + p;
         }
         return func;
     }
     ~NewtonianInterpolator() {
         for (unsigned int i = 0; i < size; i++) delete matrix[i];
         delete[]matrix;
     }
};

int main()
{
    std::vector <double> ox, error;
    double step = 0.125;
    for (int j = 0; j < 7; j++) {
        std::vector <double> x, vect;
        double l = pow(2, j) * step;
        double pi = 3.14159265358979323846;
        ox.push_back(l);
        
        int N = 5;
        for (int i = 0; i < N; i++) {
            x.push_back(i*l/(N-1));
            vect.push_back(cos(x[i]));
        }
        
        NewtonianInterpolator A(x, vect);
        
        double delta_max = 0;
        for (int k = 0; k < 1000; k++) {
            double theory = cos(k * l / 999);
            double real = A.interpolate(k * l / 999);
            double delta = theory - real;
            if (delta < 0) delta = delta * (-1); 
            if (delta > delta_max) delta_max = delta;
        } 
        error.push_back(delta_max);
        
    }
    for (int i = 0; i < 7; i++) {
        std::cout << error[i] << ' ';
    }
}