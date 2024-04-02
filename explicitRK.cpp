#include <iostream>
#include <Eigen/Dense>
#include <array>
#include <functional>
#include <Eigen/Core>
#include <vector>
#include <algorithm>

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;
using Time = double;

/*** Структура состояния ***/
struct State {
    Vec state;
    Time t;
};

/** Функция правой части дифференциального уравнения y' = f(t, y)
 *
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */
[[nodiscard]] Vec rightPart(const Time& time, const Vec& state) noexcept
{

}

/** Таблица Бутчера **/
template <unsigned s>
struct ButcherTable {
    std::array<double, s> column;
    std::array<double, s> string;
    std::array<std::array<double, s>, s> matrix;
};

/** Функция явного метода рунге-Кутты
 *
 * @tparam s  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица бутчера метода
 * @return массив решений
 */
template <unsigned s>
[[nodiscard]] std::vector<Vec> explicitRK(const State& initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)>& rightPart,
    const ButcherTable<s>& table) noexcept
{
    std::vector <Vec> K;
    Time t = initial.t;
    Vec y = initial.state;
    std::vector<Vec> res(iterations, Vec::Zero(6));

    for (unsigned i = 0; i < iterations; i++)
    {
        Vec K1 = rightPart(t, y);
        Vec K2 = rightPart(t + table.column[1] * step, y + K1 * table.matrix[1][0] * step));
        Vec K3 = rightPart(t + table.column[2] * step, y + K1 * table.matrix[2][0] * step
            + K2 * table.matrix[2][1] * step));
            Vec K4 = rightPart(t + table.column[3] * step, y + K1 * table.matrix[3][0] * step
                + K2 * table.matrix[3][1] * step
                + K3 * table.matrix[3][2] * step));

                y = y + (K1 * table.string[0] + K2 * table.string[1]
                    + K3 * table.string[2] + K4 * table.string[3]) * step;
                t = t + step;
                res[i] = y;
    }
    return(res);
}

class CSVWriter
{
    std::string fileName;
    std::string delimeter;
    int linesCount;
public:
    CSVWriter(std::string filename, std::string delm = " ") :
        fileName(filename), delimeter(delm), linesCount(0)
    {}
    template<typename T>
    void addDatainRow(T first, T last);
};

template<typename T>
void CSVWriter::addDatainRow(T first, T last)
{
    std::fstream file;
    // Open the file in truncate mode if first line else in Append Mode
    file.open(fileName, std::ios::out | (linesCount ? std::ios::app : std::ios::trunc));
    // Iterate over the range and add each lement to file seperated by delimeter.
    for (; first != last; )
    {
        file << *first;
        if (++first != last)
            file << delimeter;
    }
    file << "\n";
    linesCount++;
    // Close the file
    file.close();
}

[[nodiscard]] Vec rightPart(const Time& time, const Vec& state)
{

    double r = sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    double C = -3.986004415e14 / r / r / r;

    Vec f = Vec::Zero(6);
    f[0] = state[3];
    f[1] = state[4];
    f[2] = state[5];
    f[3] = C * state[0];
    f[4] = C * state[1];
    f[5] = C * state[2];
    return f;
}

double error(double step, ButcherTable<4> table)
{
    Vec pos{ {6378.1e3, 0., 0} };
    double v = sqrt(3.9860044158e14 / pos.norm());
    double omega = v / pos.norm();
    Vec vel{ {0, v, 0} };
    Vec y{ {pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]} };
    double max_err = 0;
    State initial;
    initial.t = 0;
    initial.state = y;
    std::vector<Vec> res = explicitRK(initial, step, 100, rightPart, table);
    double delta = step * 100;
    double y_an = pos[0] * sin(omega * delta);
    double x_an = pos[0] * cos(omega * delta);

    double error = fabs(res.back()[0] - x_an) + fabs(res.back()[1] - y_an);
    std::cout << res.back() << std::endl;
    max_err = std::max(max_err, error);
    return max_err;
}

int main() {
    const int s = 4;

    ButcherTable <s> table;

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            if (i <= j) table.matrix[i][j] = 0;
        }
    }

    table.matrix[0][0] = 0;
    table.matrix[1][0] = 0.5;
    table.matrix[1][1] = 0;
    table.matrix[2][0] = 0;
    table.matrix[2][1] = 0.5;
    table.matrix[3][0] = 0;
    table.matrix[3][1] = 0;
    table.matrix[3][2] = 1;
    table.matrix[3][3] = 0;

    table.column[0] = 0;
    table.column[1] = 0.5;
    table.column[2] = 0.5;
    table.column[3] = 1;
    table.string[0] = 1.0 / 6;
    table.string[1] = 1.0 / 3;
    table.string[2] = 1.0 / 3;
    table.string[3] = 1.0 / 6;

    Vec pos{ {6378.1e3, 0., 0} };
    double v = sqrt(3.986004415e14 / sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]));
    Vec vel{ {0, v, 0} };
    Vec y{ {pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]} };
    unsigned int step = 100;
    State initial;
    initial.t = 0;
    initial.state = y;

    Vec err = Vec::Zero(8);
    Vec steps = Vec::Zero(8);
    int i = 0;
    for (double k = 1; k <= 100; k *= 2)
    {
        i++;
        steps[i] = k;
        err[i] = error(k, table);
    }
    CSVWriter steps_csv("steps.csv");
    steps_csv.addDatainRow(steps.begin(), steps.end());
    CSVWriter mistakes_csv("mistakes.csv");
    mistakes_csv.addDatainRow(err.begin(), err.end());
    return 0;
}