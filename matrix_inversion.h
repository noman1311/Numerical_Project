
#ifndef MATRIX_INVERSION_H
#define MATRIX_INVERSION_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

#define MAX_DIM 5 // Maximum dimension of the matrix

// Function to perform matrix inversion using Gauss-Jordan elimination
void matrix_inversion() {
    int dimension;

    // User input for matrix dimension
    cout << "Enter the dimension of the square matrix (1 to 5): ";
    cin >> dimension;

    if (dimension < 1 || dimension > MAX_DIM) {
        cout << "Error: The matrix dimension must be between 1 and 5." << endl;
        return;
    }

    // Input matrix
    vector<vector<double>> matrix(dimension, vector<double>(dimension));
    vector<vector<double>> inverse(dimension, vector<double>(dimension));

    cout << "Enter the elements of the matrix row by row:" << endl;
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            cin >> matrix[i][j];
        }
    }

    // Augment the original matrix with the identity matrix of the same dimension
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (i == j)
                inverse[i][j] = 1.0;
            else
                inverse[i][j] = 0.0;
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < dimension; i++) {
        // Find the maximum element in the current column to use as pivot
        double max_element = fabs(matrix[i][i]);
        int max_row = i;

        for (int k = i + 1; k < dimension; k++) {
            if (fabs(matrix[k][i]) > max_element) {
                max_element = fabs(matrix[k][i]);
                max_row = k;
            }
        }

        // Swap the current row with the row having the maximum element
        if (max_row != i) {
            swap(matrix[i], matrix[max_row]);
            swap(inverse[i], inverse[max_row]);
        }

        // Check if the matrix is singular
        if (fabs(matrix[i][i]) < 1e-9) {
            cout << "The matrix is singular and cannot be inverted." << endl;
            return;
        }

        // Make the pivot element equal to 1
        double pivot = matrix[i][i];
        for (int j = 0; j < dimension; j++) {
            matrix[i][j] /= pivot;
            inverse[i][j] /= pivot;
        }

        // Make all elements in the current column (except the pivot) equal to 0
        for (int k = 0; k < dimension; k++) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int j = 0; j < dimension; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                    inverse[k][j] -= factor * inverse[i][j];
                }
            }
        }
    }

    // Output the inverse matrix
    cout << "The inverse of the matrix is:" << endl;
    cout << fixed << setprecision(6);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            cout << setw(10) << inverse[i][j] << " ";
        }
        cout << endl;
    }
}

#endif // MATRIX_INVERSION_H
