#include <iostream>
#define USE_MATH_DEFINES
#include <cmath>
#include <vector>

double f(double x, double y) {
    return (-0.7 * x - 0.7 * y);
}

double Runge(double t, double y0, double step) {
    double k1, k2, k3, k4;
    k1 = f(t, y0);
    k2 = f(t + step / 2, y0 + step * k1 / 2);
    k3 = f(t + step / 2, y0 + step * k2 / 2);
    k4 = f(t + step, y0 + step * k3);
    double y = y0 + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
    return y;
}

int main() {
    double initVal = 1;
    double step = 0.01;
    std::vector<double>time, y;
    time.push_back(0);
    y.push_back(initVal);
    for (double t = step; t < 1; t+=step) {
        double dy = Runge(time.back(), y.back(), step);
        time.push_back(t);
        y.push_back(dy);
    }
    double C = 1 + (-0.7 / (0.7 * 0.7));
    std::vector<double> u;
    for (double i : time) {
        double ut = -0.7 / 0.7 * (i - 1 / 0.7) + C * exp(-0.7 * i);
        u.push_back(ut);
    }
    for (int i = 0; i < time.size(); i++) {
        std::cout << "u(t) = " << u[i] << " y(t) = " << y[i] << " |u(t) - y(t)| = " << fabs(u[i] - y[i]) << std::endl;
    }
    double maxdev = 0;
    for (int i = 0; i < time.size(); i++) {
        if (maxdev < fabs(u[i] - y[i])) {
            maxdev = fabs(u[i] - y[i]);
        }
    }
    std::cout << "Max|u(t) - y(t)| = " << maxdev << std::endl;
    return 0;
}
