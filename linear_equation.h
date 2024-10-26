#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

#define MAX_DIM 5

bool diagonal_dom(vector<vector<double>>& A) {
    for (size_t i = 0; i < A.size(); i++) {
        double sum = 0.0;
        for (size_t j = 0; j < A[i].size(); j++) {
            if (j == i) continue;
            sum += fabs(A[i][j]);
        }
        if (sum >= fabs(A[i][i])) {
            return false;
        }
    }
    return true;
}

void jacobi() {
    int eqn_num;
    cout << "Enter number of equations: ";
    cin >> eqn_num;
    double tolerance;
    int maxIteration;
    cout << "Enter tolerance and max iterations: ";
    cin >> tolerance >> maxIteration;

    vector<vector<double>> A(eqn_num, vector<double>(eqn_num));
    vector<double> B(eqn_num);
    vector<double> x(eqn_num, 0.0);
    vector<double> x_new(eqn_num, 0.0);

    cout << "Enter the coefficient matrix: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        for (int j = 0; j < eqn_num; j++) {
            cin >> A[i][j];
        }
    }

    if (!diagonal_dom(A)) {
        cout << "The matrix is not diagonally dominant." << endl;
        return;
    }

    cout << "Enter the constant matrix: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        cin >> B[i];
    }

    double MAE = 0.0;
    for (int iter = 0; iter < maxIteration; iter++) {
        for (int i = 0; i < eqn_num; i++) {
            double sum = 0;
            for (int j = 0; j < eqn_num; j++) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (B[i] - sum) / A[i][i];
        }

        double error = 0.0;
        for (int i = 0; i < eqn_num; i++) {
            error += fabs(x_new[i] - x[i]);
        }
        MAE = error / eqn_num;

        if (MAE <= tolerance) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
        x = x_new;
    }

    if (MAE > tolerance) {
        cout << "The roots did not converge with the given number of iterations." << endl;
        return;
    }

    for (int i = 0; i < eqn_num; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

void gaussSeidel() {
    int eqn_num;
    cout << "Enter number of equations: ";
    cin >> eqn_num;
    int maxIteration;
    double tolerance;
    cout << "Enter max iterations and tolerance: ";
    cin >> maxIteration >> tolerance;

    vector<vector<double>> A(eqn_num, vector<double>(eqn_num));
    vector<double> B(eqn_num);
    vector<double> x(eqn_num, 0.0);

    cout << "Enter the coefficient matrix: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        for (int j = 0; j < eqn_num; j++) {
            cin >> A[i][j];
        }
    }

    if (!diagonal_dom(A)) {
        cout << "The matrix is not diagonally dominant." << endl;
        return;
    }

    cout << "Enter the constant matrix: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        cin >> B[i];
    }

    for (int iter = 0; iter < maxIteration; iter++) {
        double error = 0.0;
        for (int i = 0; i < eqn_num; i++) {
            double sum = 0;
            for (int j = 0; j < eqn_num; j++) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            double x_new = (B[i] - sum) / A[i][i];
            error += fabs(x_new - x[i]);
            x[i] = x_new;
        }

        double MAE = error / eqn_num;
        if (MAE <= tolerance) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            for (int i = 0; i < eqn_num; i++) {
                cout << "x" << i + 1 << " = " << x[i] << endl;
            }
            return;
        }
    }
    cout << "Did not converge within the maximum number of iterations." << endl;
}

void gauss_elimination() {
    int eqn_num;
    cout << "Enter number of equations: ";
    cin >> eqn_num;

    vector<vector<double>> A(eqn_num, vector<double>(eqn_num + 1));
    cout << "Enter the augmented matrix: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        for (int j = 0; j <= eqn_num; j++) {
            cin >> A[i][j];
        }
    }

    // Forward elimination
    for (int i = 0; i < eqn_num - 1; i++) {
        for (int k = i + 1; k < eqn_num; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= eqn_num; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Back substitution
    vector<double> x(eqn_num);
    for (int i = eqn_num - 1; i >= 0; i--) {
        x[i] = A[i][eqn_num];
        for (int j = i + 1; j < eqn_num; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    cout << "Solution: " << endl;
    for (int i = 0; i < eqn_num; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

void gauss_jordan() {
    int num_variables;
    cout << "Enter the number of variables (1 to 5): ";
    cin >> num_variables;

    if (num_variables < 1 || num_variables > MAX_DIM) {
        cout << "Error: Number of variables must be between 1 and 5." << endl;
        return;
    }

    vector<vector<double>> augmented_matrix(num_variables, vector<double>(num_variables + 1, 0));
    cout << "Enter the coefficients and constants for each equation (augmented matrix format):" << endl;
    for (int i = 0; i < num_variables; i++) {
        for (int j = 0; j <= num_variables; j++) {
            cin >> augmented_matrix[i][j];
        }
    }

    for (int i = 0; i < num_variables; i++) {
        double max_element = fabs(augmented_matrix[i][i]);
        int max_row = i;

        for (int k = i + 1; k < num_variables; k++) {
            if (fabs(augmented_matrix[k][i]) > max_element) {
                max_element = fabs(augmented_matrix[k][i]);
                max_row = k;
            }
        }

        if (max_row != i) {
            swap(augmented_matrix[i], augmented_matrix[max_row]);
        }

        if (fabs(augmented_matrix[i][i]) < 1e-9) {
            cout << "The system has no unique solution (singular matrix)." << endl;
            return;
        }

        double pivot = augmented_matrix[i][i];
        for (int j = 0; j <= num_variables; j++) {
            augmented_matrix[i][j] /= pivot;
        }

        for (int k = 0; k < num_variables; k++) {
            if (k != i) {
                double factor = augmented_matrix[k][i];
                for (int j = 0; j <= num_variables; j++) {
                    augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
                }
            }
        }
    }

    cout << "Solution:" << endl;
    cout << fixed << setprecision(6);
    for (int i = 0; i < num_variables; i++) {
        cout << "x" << i + 1 << " = " << augmented_matrix[i][num_variables] << endl;
    }
}

bool perform_lu_factorization(vector<vector<double>>& matrix, vector<vector<double>>& lower_matrix, vector<vector<double>>& upper_matrix) {
    int n = matrix.size();
    lower_matrix = vector<vector<double>>(n, vector<double>(n, 0));
    upper_matrix = vector<vector<double>>(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += lower_matrix[i][j] * upper_matrix[j][k];
            }
            upper_matrix[i][k] = matrix[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k) {
                lower_matrix[i][i] = 1.0;
            } else {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += lower_matrix[k][j] * upper_matrix[j][i];
                }
                if (upper_matrix[i][i] == 0) {
                    return false;
                }
                lower_matrix[k][i] = (matrix[k][i] - sum) / upper_matrix[i][i];
            }
        }
    }
    return true;
}


void solve_lu() {
    int num_variables;
    cout << "Enter the number of variables (1 to 5): ";
    cin >> num_variables;

    if (num_variables < 1 || num_variables > MAX_DIM) {
        cout << "Error: Number of variables must be between 1 and 5." << endl;
        return;
    }

    vector<vector<double>> augmented_matrix(num_variables, vector<double>(num_variables + 1, 0));
    cout << "Enter the coefficients and constants for each equation (augmented matrix format):" << endl;
    for (int i = 0; i < num_variables; i++) {
        for (int j = 0; j <= num_variables; j++) {
            cin >> augmented_matrix[i][j];
        }
    }

    vector<vector<double>> lower_matrix, upper_matrix;
    vector<vector<double>> coefficient_matrix(num_variables, vector<double>(num_variables, 0));
    vector<double> constants(num_variables, 0);

    for (int i = 0; i < num_variables; i++) {
        for (int j = 0; j < num_variables; j++) {
            coefficient_matrix[i][j] = augmented_matrix[i][j];
        }
        constants[i] = augmented_matrix[i][num_variables];
    }

    if (!perform_lu_factorization(coefficient_matrix, lower_matrix, upper_matrix)) {
        cout << "Error: LU Factorization failed." << endl;
        return;
    }

    vector<double> y(num_variables, 0);
    for (int i = 0; i < num_variables; i++) {
        y[i] = constants[i];
        for (int j = 0; j < i; j++) {
            y[i] -= lower_matrix[i][j] * y[j];
        }
    }

    vector<double> x(num_variables, 0);
    for (int i = num_variables - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < num_variables; j++) {
            x[i] -= upper_matrix[i][j] * x[j];
        }
        x[i] /= upper_matrix[i][i];
    }

    cout << "Solution:" << endl;
    cout << fixed << setprecision(6);
    for (int i = 0; i < num_variables; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}



void console_linear_equation() {
    while(true){
        int choice;
        cout << "Select a method to solve the system of linear equations:" << endl;
        cout << "1. Jacobi Iterative Method" << endl;
        cout << "2. Gauss-Seidel Iterative Method" << endl;
        cout << "3. Gauss Elimination Method" << endl;
        cout << "4. Gauss-Jordan Elimination Method" << endl;
        cout << "5. LU Factorization Method" << endl;
        cout << "Enter your choice (1-5): ";
        cin >> choice;

        switch (choice) {
            case 1:
                jacobi();
                break;
            case 2:
                gaussSeidel();
                break;
            case 3:
                gauss_elimination();
                break;
            case 4:
                gauss_jordan();
                break;
            case 5:
                solve_lu();
                break;
            default:
                cout << "Invalid choice. Please enter a number between 1 and 5.\n\n" << endl;
        }
        cout<<"Enter :\n 1 to stay in this menu\n 0 to go back to main menu\n";
        int ttt;
        cin>>ttt;
        if(ttt==0)
            break;
    }
}

#endif // LINEAR_SYSTEM_SOLVER_H
