#pragma once
void mygen(double **a, double **a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema);
void Q_matrix(double** Q, int n, int schema);
void matr_mul(double **a, double **b, double **c, int n);
double matr_inf_norm(double **a, int n);
double v_inf_norm(double *v, int n);
//*///*//

void Print(double** a, int n);

// Ваши функции
int your_lu_inv(double **a, int n);     // обращение матрицы
int your_lu_solv(double **a, int n, double *b);  // решение слау

