#include "progr.h"
#include <cmath>
#include <iomanip>
#include <iostream>


using namespace std;

double **createIdentityMatrix(int n) {
  double **e = new double *[n];
  for (int i = 0; i < n; ++i) {
    e[i] = new double[n];
  }
  for (int i = 0; i < n; ++i) {
    e[i][i] = 1;
  }
  return e;
}

double **inverseMatrix(double **matrix, int n) {
  // Создаем единичную матрицу
  double **identity = createIdentityMatrix(n);

  for (int i = 0; i < n; i++) {
    // Поиск максимального элемента в столбце
    int maxRow = i;
    for (int j = i + 1; j < n; j++) {
      if (abs(matrix[j][i]) > abs(matrix[maxRow][i])) {
        maxRow = j;
      }
    }

    // Перестановка строк, если максимальный элемент не на диагонали
    if (maxRow != i) {
      swap(matrix[i], matrix[maxRow]);
      swap(identity[i], identity[maxRow]);
    }

    // Проверка того что главный элемент не равен 0, иначе определитель 0
    double pivot = matrix[i][i];
    // if (pivot < 0.000001) {
    //   cout << "Определитель матрицы равен 0!!!\n";
    //   exit(1);
    // }

    // Деление строки на главный элемент
    for (int j = i; j < n; j++) {
      matrix[i][j] /= pivot;
    }
    for (int j = 0; j < n; j++) {
      identity[i][j] /= pivot;
    }

    // Вычитание строк
    for (int k = 0; k < n; k++) {
      if (k != i) {
        double factor = matrix[k][i];
        for (int j = i; j < n; j++) {
          matrix[k][j] -= factor * matrix[i][j];
        }
        for (int j = 0; j < n; j++) {
          identity[k][j] -= factor * identity[i][j];
        }
      }
    }
  }

  return identity;
}

int main() {
  int n;
  double **matrix;
  double **matrix_copy;
  double **inverse1;
  double **inc;
  double normA;
  double normA_inv;
  double obusl;
  double normR_gen;
  double normR;
  double alpha;
  double beta;
  cout << "|" << setw(6) << "alpha" << "|" << setw(6) << "beta" << "|"
       << setw(11) << " ||  A  || " << "|" << setw(11) << " ||A_inv|| " << "|"
       << setw(11) << "   obus    " << "|" << setw(11) << " ||  Z  || " << "|"
       << setw(11) << "   dzeta   " << "|" << setw(11) << " ||  r  || " << "|"
       << setw(11) << "   R_gen   " << "|\n"
       << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
          "|||||||||||||||||||||||||||||||\n";

  alpha = 1;
  beta = 1;

  for (int j = 1; j <= 18; j++) {

    n = 100;
    matrix = new double *[n];
    for (int i = 0; i < n; ++i) {
      matrix[i] = new double[n];
    }

    inverse1 = new double *[n];
    for (int i = 0; i < n; ++i) {
      inverse1[i] = new double[n];
    }

    mygen(matrix, inverse1, n, alpha, beta, 1, 2, 0, 1, &normA, &normA_inv,
          &obusl, &normR_gen);

    // cout << matrix[0][0] << endl;
    // beta = pow(10, j);
    alpha /= 10;
    matrix_copy = new double *[n];
    for (int i = 0; i < n; ++i) {
      matrix_copy[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
      for (int l = 0; l < n; l++) {
        matrix_copy[i][l] = matrix[i][l];
      }
    }

    // Print(matrix, n);
    //  Нахождение обратной матрицы
    inverseMatrix(matrix, n);
    // Print(matrix, n);

    inc = new double *[n];
    for (int i = 0; i < n; ++i) {
      inc[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
      for (int l = 0; l < n; l++) {
        inc[i][l] = inverse1[i][l] - matrix[i][l];
      }
    }
    normR = matr_inf_norm(inc, n);
    /*cout << "||R|| = " << normR << endl;*/

    double **r = new double *[n];
    for (int i = 0; i < n; i++)
      r[i] = new double[n];
    matr_mul(matrix_copy, matrix, r, n);
    for (int i = 0; i < n; i++)
      r[i][i] -= 1.;
    double myNormR = matr_inf_norm(r, n);

    /*printMatrix(r, n);
    cout << "||||\n";*/
    // printMatrix(inverse1, n);
    /*cout << "||||\n";
    printMatrix(inverse, n);*/
    cout << "|" << setw(6) << alpha << "|" << setw(6) << beta << "|" << setw(11)
         << normA << "|" << setw(11) << normA_inv << "|" << setw(11) << obusl
         << "|" << setw(11) << normR << "|" << setw(11) << normR / normA_inv
         << "|" << setw(11) << myNormR << "|" << setw(11) << normR_gen << "|\n"
         << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
            "|||||||||||||||||||||||||||||||||\n";

    // Освобождение памяти
    for (int i = 0; i < n; ++i) {
      delete[] matrix[i];
      delete[] matrix_copy[i];
      delete[] inverse1[i];
      delete[] inc[i];
    }
    delete[] matrix;
    delete[] matrix_copy;
    delete[] inverse1;
    delete[] inc;
  }
  return 0;
}
