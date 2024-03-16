#include <iostream>
#include <vector>

using namespace std;

const double EPSILON = 1e-10;

// Perform Gaussian elimination
void gaussianElimination(vector<vector<double>> &A) {
  int n = A.size();

  for (int i = 0; i < n; i++) {
    // Find the pivot row
    int pivotRow = i;
    for (int j = i + 1; j < n; j++) {
      if (abs(A[j][i]) > abs(A[pivotRow][i]))
        pivotRow = j;
    }

    // Swap rows
    if (pivotRow != i) {
      swap(A[i], A[pivotRow]);
    }

    // Make the diagonal elements 1
    double divisor = A[i][i];
    for (int j = 0; j < n; j++) {
      A[i][j] /= divisor;
    }

    // Eliminate non-zero elements below the diagonal
    for (int j = i + 1; j < n; j++) {
      double factor = A[j][i];
      for (int k = 0; k < n; k++) {
        A[j][k] -= factor * A[i][k];
      }
    }
  }
}

// Perform back substitution to find the inverse
void backSubstitution(vector<vector<double>> &A) {
  int n = A.size();

  for (int i = n - 1; i >= 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      double factor = A[j][i];
      for (int k = 0; k < n; k++) {
        A[j][k] -= factor * A[i][k];
      }
    }
  }
}

// Function to find inverse matrix
void inverseMatrix(vector<vector<double>> &matrix) {
  int n = matrix.size();

  // Perform Gaussian elimination
  gaussianElimination(matrix);

  // Perform back substitution
  backSubstitution(matrix);
}

// Function to print matrix
void printMatrix(const vector<vector<double>> &matrix) {
  for (const auto &row : matrix) {
    for (double element : row) {
      cout << element << " ";
    }
    cout << endl;
  }
}

int main() {
  vector<vector<double>> matrix = {{1, 3, 3}, {1, 4, 3}, {1, 3, 4}};
  cout << "Original Matrix:" << endl;
  printMatrix(matrix);

  inverseMatrix(matrix);
  cout << "Inverse Matrix:" << endl;
  printMatrix(matrix);

  return 0;
}
