#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Define the differential equation dy/dx = c1 * x^n1 + c2 * y^n2
double differential_equation(double x, double y, double c1, double n1, double c2, double n2) {
    return c1 * pow(x, n1) + c2 * pow(y, n2);
}

// Runge-Kutta 4th Order Method function
void runge_kutta() {
    double initial_x, initial_y, target_x, step_size;
    double coefficient1, power1, coefficient2, power2;

    // Input initial conditions and step size
    cout << "Enter the initial value of x (x0): ";
    cin >> initial_x;
    cout << "Enter the initial value of y (y0): ";
    cin >> initial_y;
    cout << "Enter the value of x at which you want to estimate y: ";
    cin >> target_x;
    cout << "Enter the step size: ";
    cin >> step_size;

    // Input the coefficients and powers for the differential equation
    cout << "Enter the coefficient c1 (for x^n1 term): ";
    cin >> coefficient1;
    cout << "Enter the power n1 (for x^n1 term): ";
    cin >> power1;
    cout << "Enter the coefficient c2 (for y^n2 term): ";
    cin >> coefficient2;
    cout << "Enter the power n2 (for y^n2 term): ";
    cin >> power2;

    // Number of iterations
    int iterations = static_cast<int>((target_x - initial_x) / step_size);
    double k1, k2, k3, k4;
    double y = initial_y;

    // Runge-Kutta 4th Order Method
    for (int i = 0; i < iterations; i++) {
        k1 = step_size * differential_equation(initial_x, y, coefficient1, power1, coefficient2, power2);
        k2 = step_size * differential_equation(initial_x + step_size / 2.0, y + k1 / 2.0, coefficient1, power1, coefficient2, power2);
        k3 = step_size * differential_equation(initial_x + step_size / 2.0, y + k2 / 2.0, coefficient1, power1, coefficient2, power2);
        k4 = step_size * differential_equation(initial_x + step_size, y + k3, coefficient1, power1, coefficient2, power2);

        y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        initial_x += step_size;
    }

    // Output the result
    cout << fixed << setprecision(6);
    cout << "The estimated value of y at x = " << target_x << " is: " << y << endl;
}

#endif // RUNGE_KUTTA_H

