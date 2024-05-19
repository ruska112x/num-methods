/*?????-?? ??? ?? main - ?????? ????????????? ?????????? ??????


#include "progr.h"
#include <iostream>
#include <math.h>
using namespace std;

#define N 10
#define ALPHA 1.
#define BETA 1e1

// ---------------------------------------------------

//ofstream logfile("log.txt"); //???? ??? ???????


int main()
{
	int i, j;

	int n = N;
	double** a = new double* [n];
	for (i = 0; i < n; i++)
		a[i] = new double[n];

	double** a_inv = new double* [n];
	for (i = 0; i < n; i++)
		a_inv[i] = new double[n];

	double alpha = ALPHA;
	double beta = BETA;

	
	mygen(a, a_inv, n, alpha, beta, 1, 2, 1, 1); //??????? ?????????
	//mygen(a, a_inv, n, alpha, beta, 0, 0, 2, 1); // ????????? ?????? (n>2)
	//mygen(a, a_inv, n, alpha, beta, 1, 2, 0, 1); // ????????????

	*/
//-------------------------------------------------------------------------


#include "progr/progr.h"
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

//ofstream logfile("log.txt"); //???? ??? ???????
//extern ofstream logfile;

ofstream logfile("log.txt"); //???? ??? ???????
							 //extern ofstream logfile;

int your_lu_inv(double** a, int n) 
{
	// ????????? ???????
	// ??? ???	
	return 1;
}


int your_lu_solv(double **a, int n, double *b)
{
	// ??????? ????
	// ??? ???	
	return 1;
}


void mygen(double **a, double **a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema, double *normA, double *normA_inv, double *obusl, double *normR_gen)
{
	int i, j, k;


	//cout << "   M A T R I X  G E N.  " << endl;

	//cout << "              N = " << n << endl;
	//cout << " | lambda_min | = " << alpha << endl;
	//cout << " | lambda_max | = " << beta << endl;

	double *lambda = new double[n];

	// ????????????? ??????
	//cout << " sign_law = " << sign_law << endl;

	double *sign = new double[n];
	for (i = 0; i<n; i++) sign[i] = 1.;

	switch (sign_law)
	{
	case -1:
		for (i = 0; i<n; i++) sign[i] = -1.;
		break;

	case 0:
		sign[0] = 1.;
		for (i = 1; i<n; i++) sign[i] = -sign[i - 1];
		break;

		//?????? ?????? ????????????? ??????
		// ...

	}
	/*	for( i=0; i<n; i++ ) cout<<sign[i]<<" ";
	cout<<endl;*/

	//????????????? ???????????? ?????
	//cout << " lambda_law = " << lambda_law << endl;

	double *kappa = new double[n];
	for (i = 0; i<n; i++) kappa[i] = (double)i / double(n - 1);
	switch (lambda_law)
	{
	case 1:
		cout << " kappa = sqrt( ) " << endl;
		for (i = 0; i<n; i++) kappa[i] = sqrt(kappa[i]);
		break;

	case 2:
		//cout << " kappa = sin( ) " << endl;
		double pi_half = acos(-1.)*0.5;
		for (i = 0; i<n; i++) kappa[i] = sin(pi_half*kappa[i]);
		break;

		//?????? ?????? ????????????? ??????????? ?????
		// ...

	}
	/*	for( i=0; i<n; i++ ) cout<<kappa[i]<<" ";
	cout<<endl;
	*/


	double *J = new double[n];
	for (i = 0; i<n; i++) J[i] = sign[i] * ((1. - kappa[i])*alpha + kappa[i] * beta);
	/*	for( i=0; i<n; i++ ) cout<<J[i]<<" ";
	cout<<endl;
	*/

	double *J_inv = new double[n];
	for (i = 0; i<n; i++) J_inv[i] = 1. / J[i];

	double **Q = new double *[n];
	for (i = 0; i<n; i++) Q[i] = new double[n];

	double aa[3];


	//cout << " variant = " << variant << endl;

	switch (variant)
	{
	case 0: //???????????? ???????
		//cout << " simmetric matrix:" << endl;
		//cout << " schema = " << schema << endl;
		switch (schema)
		{
		case 1:
			Q_matrix(Q, n, schema);

			for (a[0][0] = 0., k = 0; k<n; k++) a[0][0] += Q[0][k] * J[k] * Q[0][k];
			for (j = 1; j<n; j++)
			{
				for (a[0][j] = 0., k = j - 1; k<n; k++) a[0][j] += Q[0][k] * J[k] * Q[j][k];
				a[j][0] = a[0][j];
			}
			for (i = 1; i<n; i++)
			{
				for (a[i][i] = 0., k = i - 1; k<n; k++) a[i][i] += Q[i][k] * J[k] * Q[i][k];
				for (j = i + 1; j<n; j++)
				{
					for (a[i][j] = 0., k = j - 1; k<n; k++) a[i][j] += Q[i][k] * J[k] * Q[j][k];
					a[j][i] = a[i][j];
				}
			}

			//_______
			for (a_inv[0][0] = 0., k = 0; k<n; k++) a_inv[0][0] += Q[0][k] * J_inv[k] * Q[0][k];
			for (j = 1; j<n; j++)
			{
				for (a_inv[0][j] = 0., k = j - 1; k<n; k++) a_inv[0][j] += Q[0][k] * J_inv[k] * Q[j][k];
				a_inv[j][0] = a_inv[0][j];
			}
			for (i = 1; i<n; i++)
			{
				for (a_inv[i][i] = 0., k = i - 1; k<n; k++) a_inv[i][i] += Q[i][k] * J_inv[k] * Q[i][k];
				for (j = i + 1; j<n; j++)
				{
					for (a_inv[i][j] = 0., k = j - 1; k<n; k++) a_inv[i][j] += Q[i][k] * J_inv[k] * Q[j][k];
					a_inv[j][i] = a_inv[i][j];
				}
			}
			break;

		}//schema
		break;

	case 1: //??????? ??????? ?????????
		cout << " simple structure matrix:" << endl;
		cout << " schema = " << schema << endl;
		switch (schema)
		{
		case 1:
			//TJ
			//?????? ??????
			a[0][0] = J[0];
			a[0][1] = -J[1];
			//			for(j=2; j<n; j++ ) a[0][j] = 0.;
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				//				for(j=0; j<i-1; j++ ) a[i][j] = 0.;  
				a[i][i - 1] = -J[i - 1];
				a[i][i] = J[i] + J[i];
				a[i][i + 1] = -J[i + 1];
				//				for(j=i+2; j<n; j++ ) a[i][j] = 0.;
			}
			//????????? (n-1)
			//			for(j=0; j<n-2; j++ ) a[n-1][j] = 0.;  
			a[n - 1][n - 2] = -J[n - 2];
			a[n - 1][n - 1] = J[n - 1] + J[n - 1];

			//(TJ)T^{-1}
			//?????? ??????
			aa[1] = a[0][0];  aa[2] = a[0][1];
			a[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
			double s = aa[1] + aa[2];
			for (j = 1; j<n; j++) a[0][j] = s * (double)(n - j);
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				aa[0] = a[i][i - 1];  aa[1] = a[i][i];  aa[2] = a[i][i + 1];
				for (j = 0; j<i; j++) a[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s = aa[0] + aa[1];
				a[i][i] = s * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s += aa[2];
				for (j = i + 1; j<n; j++) a[i][j] = s * (double)(n - j);
			}
			//????????? (n-1)
			aa[0] = a[n - 1][n - 2];  aa[1] = a[n - 1][n - 1];
			s = aa[0] + aa[0] + aa[1];
			for (j = 0; j<n - 1; j++) a[n - 1][j] = s;
			a[n - 1][n - 1] = aa[0] + aa[1];
			//_______

			//TJ^{-1}
			//?????? ??????
			a_inv[0][0] = J_inv[0];
			a_inv[0][1] = -J_inv[1];
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				a_inv[i][i - 1] = -J_inv[i - 1];
				a_inv[i][i] = J_inv[i] + J_inv[i];
				a_inv[i][i + 1] = -J_inv[i + 1];
			}
			//????????? (n-1)
			a_inv[n - 1][n - 2] = -J_inv[n - 2];
			a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];

			//(TJ^{-1})T^{-1}
			//?????? ??????
			aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
			a_inv[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
			s = aa[1] + aa[2];
			for (j = 1; j<n; j++) a_inv[0][j] = s * (double)(n - j);
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				aa[0] = a_inv[i][i - 1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i + 1];
				for (j = 0; j<i; j++) a_inv[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s = aa[0] + aa[1];
				a_inv[i][i] = s * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s += aa[2];
				for (j = i + 1; j<n; j++) a_inv[i][j] = s * (double)(n - j);
			}
			//????????? (n-1)
			aa[0] = a_inv[n - 1][n - 2];  aa[1] = a_inv[n - 1][n - 1];
			s = aa[0] + aa[0] + aa[1];
			for (j = 0; j<n - 1; j++) a_inv[n - 1][j] = s;
			a_inv[n - 1][n - 1] = aa[0] + aa[1];
			break;

		}//schema
		break;

	case 2: //???? ????????? ?????? 2x2 ??? ??????????? ?.?.
		cout << " J_2 type matrix:" << endl;
		cout << " schema = " << schema << endl;

		switch (schema)
		{
		case 1:
			//TJ
			//?????? ??????
			a[0][0] = J[0];
			a[0][1] = 1. - J[0];
			//?????? ??????
			a[1][0] = -J[0];
			a[1][1] = -1. + J[0] + J[0];
			a[1][2] = -J[2];
			//?????? ??????
			a[2][1] = -J[0];
			a[2][2] = J[2] + J[2];
			if (n>3) a[2][3] = -J[3];
			//?? ?????????
			for (i = 3; i<n - 1; i++)
			{
				a[i][i - 1] = -J[i - 1];
				a[i][i] = J[i] + J[i];
				a[i][i + 1] = -J[i + 1];
			}
			//????????? (n-1)
			if (n>3)
			{
				a[n - 1][n - 2] = -J[n - 2];
				a[n - 1][n - 1] = J[n - 1] + J[n - 1];
			}

			//(TJ)T^{-1}
			//?????? ??????
			aa[1] = a[0][0];  aa[2] = a[0][1];
			a[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
			double s = aa[1] + aa[2];
			for (j = 1; j<n; j++) a[0][j] = s * (double)(n - j);
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				aa[0] = a[i][i - 1];  aa[1] = a[i][i];  aa[2] = a[i][i + 1];
				for (j = 0; j<i; j++) a[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s = aa[0] + aa[1];
				a[i][i] = s * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s += aa[2];
				for (j = i + 1; j<n; j++) a[i][j] = s * (double)(n - j);
			}
			//????????? (n-1)
			aa[0] = a[n - 1][n - 2];  aa[1] = a[n - 1][n - 1];
			s = aa[0] + aa[0] + aa[1];
			for (j = 0; j<n - 1; j++) a[n - 1][j] = s;
			a[n - 1][n - 1] = aa[0] + aa[1];
			//_______
			//TJ^{-1}
			//?????? ??????
			a_inv[0][0] = J_inv[0];
			a_inv[0][1] = -J_inv[0] * J_inv[0] - J_inv[0];
			//?????? ??????
			a_inv[1][0] = -J_inv[0];
			a_inv[1][1] = J_inv[0] * J_inv[0] + J_inv[0] + J_inv[0];
			a_inv[1][2] = -J_inv[2];
			//?????? ??????
			a_inv[2][1] = -J_inv[0];
			a_inv[2][2] = J_inv[2] + J_inv[2];
			if (n>3) a_inv[2][3] = -J_inv[3];
			//?? ?????????
			for (i = 3; i<n - 1; i++)
			{
				a_inv[i][i - 1] = -J_inv[i - 1];
				a_inv[i][i] = J_inv[i] + J_inv[i];
				a_inv[i][i + 1] = -J_inv[i + 1];
			}
			//????????? (n-1)
			if (n>3)
			{
				a_inv[n - 1][n - 2] = -J_inv[n - 2];
				a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];
			}

			//(TJ^{-1})T^{-1}
			//?????? ??????
			aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
			a_inv[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
			s = aa[1] + aa[2];
			for (j = 1; j<n; j++) a_inv[0][j] = s * (double)(n - j);
			//?? ?????????
			for (i = 1; i<n - 1; i++)
			{
				aa[0] = a_inv[i][i - 1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i + 1];
				for (j = 0; j<i; j++) a_inv[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s = aa[0] + aa[1];
				a_inv[i][i] = s * (double)(n - i) + aa[2] * (double)(n - i - 1);
				s += aa[2];
				for (j = i + 1; j<n; j++) a_inv[i][j] = s * (double)(n - j);
			}
			//????????? (n-1)
			aa[0] = a_inv[n - 1][n - 2];  aa[1] = a_inv[n - 1][n - 1];
			s = aa[0] + aa[0] + aa[1];
			for (j = 0; j<n - 1; j++) a_inv[n - 1][j] = s;
			a_inv[n - 1][n - 1] = aa[0] + aa[1];


			break;
		}//schema

		break;

	}//variant

	 //______________________________________________________________________


	*normA = matr_inf_norm(a, n);
	//cout << " ||  A  || = " << *norm << endl;

	*normA_inv = matr_inf_norm(a_inv, n);
	//cout << " ||A_inv|| = " << *normA_inv << endl;

	*obusl = *normA * *normA_inv;
	//cout << " obusl = " << *obusl << endl;

	//??????? ?????????
	double **r = new double*[n];
	for (i = 0; i < n; i++)
		r[i] = new double[n];
	matr_mul(a, a_inv, r, n);
	for (i = 0; i < n; i++) r[i][i] -= 1.;

	//cout << "a:" << endl;
	if (n < 6) {

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++) cout << a_inv[i][j] << "\t";
			cout << endl;
		}
	}

	*normR_gen = matr_inf_norm(r, n);
	//cout << " ||R_gen|| = " << *normR_gen << endl;


}//mygen

void Q_matrix(double **Q, int n, int schema)
{
	int i, j;
	double  q;

	double curr, next = 1.;
	for (j = 0; j<n - 1; j++)
	{
		curr = next;
		next += 1.;

		q = 1. / sqrt(curr*next);
		for (i = 0; i <= j; i++) Q[i][j] = q;
		Q[j + 1][j] = -sqrt(curr / next);
		for (i = j + 2; i<n; i++) Q[i][j] = 0.;
	}

	q = 1. / sqrt((double)n);
	for (i = 0; i<n; i++) Q[i][n - 1] = q;
}


void matr_mul(double **a, double **b, double **c, int n)
{
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (c[i][j] = 0., k = 0; k<n; k++) c[i][j] += a[i][k] * b[k][j];

}

double matr_inf_norm(double **a, int n)
{
	int i, j;
	double s, norm = 0.;

	for (i = 0; i < n; i++)
	{
		for (s = 0., j = 0; j < n; j++) s += fabs(a[i][j]);
		if (s > norm) norm = s;
	}

	return norm;
}
double v_inf_norm(double* v, int n)
{
	int i;
	double t, norm = 0.;
	for (i = 0; i < n; i++)
		if ((t = fabs(v[i])) > norm) norm = t;

	return norm;
}
void Print(double** a, int n)
{
	cout << "a:" << endl;
	if (n < 7) {

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)  cout <<  a[i][j] << "     "; 
			cout << endl;
		}
	}
}