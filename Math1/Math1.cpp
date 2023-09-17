#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

double f(double* y, double t, int i)
{
	double w;

	switch (i)
	{
	case 0: return w = 2 * (y[0] - y[0] * y[1]); break;
	case 1: return w = -(y[1] - y[0] * y[1]); break;
	}
}

const int n = 2;
const double tMax = 7.0;
double y[n] = { 1.0, 3.0 };

void EulerMethod(double (*f)(double*, double, int), int n, double tMax, double* yc)
{
	setlocale(LC_ALL, "Russian");

	const double t0 = 0.0;
	const double tau = 0.01;
	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double [n]{ 0.0 };
	double tn, tk, deltat;

	tn = omp_get_wtime();

	for (double t = t0; t <= tMax; t += tau)
	{
		for (int i = 0; i < n; ++i)
		{
			yy[i] = y[i] + tau * f(y, t, i);
		}

		for (int i = 0; i < n; ++i)
		{
			y[i] = yy[i];
		}
	}

	tk = omp_get_wtime();
	deltat = tk - tn;
	std::cout << "Метод Эйлера" << std::endl;
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{

		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}

void RK2(double (*f)(double*, double, int), int n, double tMax, double yc[])
{
	setlocale(LC_ALL, "Russian");

	const double t0 = 0.0;
	const double tau = 0.01;
	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double [n] { 0.0 };
	double tn, tk, deltat;
	double* ff = new double [n]{ 0.0 };

	tn = omp_get_wtime();
	for (double t = t0; t <= tMax; t += tau)
	{
		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + 0.5 * tau * f(y, t, i);
		}

		for (int i = 0; i < n; i++)
		{
			ff[i] = f(yy, t + 0.5 * tau, i);
		}

		for (int i = 0; i < n; i++)
		{
			y[i] += tau * ff[i];
		}
	}

	tk = omp_get_wtime();
	deltat = tk - tn;

	std::cout << "Метод Рунге-Кутты 2" << std::endl;
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{

		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}

void Correction(double (*f)(double*, double, int), int n, double tMax, double yc[])
{
	setlocale(LC_ALL, "Russian");

	const double t0 = 0.0;
	const double tau = 0.01;
	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double [n] { 0.0 };
	double tn, tk, deltat;
	double* ff = new double [n] { 0.0 };

	tn = omp_get_wtime();
	for (double t = t0; t <= tMax; t += tau)
	{
		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + tau * f(y, t, i);
		}

		for (int i = 0; i < n; i++)
		{
			ff[i] = y[i] + tau * (f(yy, t + tau, i) + f(y, t, i)) / 2.0;
		}

		for (int i = 0; i < n; i++)
		{
			y[i] = ff[i];
		}
	}

	tk = omp_get_wtime();
	deltat = tk - tn;

	std::cout << "Метод \"Прогноз коррекции\"" << std::endl;
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{
		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}

void RK4(double (*f)(double*, double, int), int n, double tMax, double yc[])
{
	setlocale(LC_ALL, "Russian");

	const int m = 4;
	const double t0 = 0.0;
	const double tau = 0.01;
	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double [n] { 0.0 };
	double tn, tk, deltat;
	double* ff = new double [n] { 0.0 };
	double** R = new double*[m]{};
	for (int i = 0; i < m; i++)
	{
		R[i] = new double[n] {0.0};
	}

	tn = omp_get_wtime();
	for (double t = t0; t < tMax; t += tau)
	{
		for (int i = 0; i < n; i++)
		{
			R[0][i] = tau * f(y, t, i);
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + 0.5 * R[0][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[1][i] = tau * f(yy, t + 0.5 * tau, i);
			
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + 0.5 * R[1][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[2][i] = tau * f(yy, t + 0.5 * tau, i);
			yy[i] = y[i] + R[2][i];
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + R[2][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[3][i] = tau * f(yy, t + tau, i);
		}

		for (int i = 0; i < n; i++)
		{
			y[i] += (R[0][i] + 2 * R[1][i] + 2 * R[2][i] + R[3][i]) / 6.0;
		}
	}

	tk = omp_get_wtime();
	deltat = tk - tn;

	std::cout << "Метод Рунге-Кутты 4" << std::endl;
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{
		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}

void ImplictEulerMethod(double (*f)(double*, double, int), int n, double tMax, double yc[])
{
	setlocale(LC_ALL, "Russian");

	const int m = 4;
	const double t0 = 0.0;
	const double tau = 0.01;
	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double[n] {0.0};
	double tn, tk, deltat;
	double deltah = 0.001, opred; 
	double* opredn = new double[n]{ 0.0 };
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n] {0.0};
	}
	double* b = new double[n] {0.0}; 
	double* p = new double[n] {0.0};

	//const int n = 2;
	//const int m = 4;
	//const double t0 = 0.0;
	//const double tMax = 7.0;
	//const double tau = 0.01;
	//double t = t0;
	//double y[n] = { 1.0, 3.0 };
	//double yy[n] = { 0.0 };
	//double tn, tk, deltat;
	//double deltah = 0.001, opred, opredn[n] = { 0.0 };
	//double a[n][n] = { 0.0 }, b[n] = { 0.0 }, p[n] = { 0.0 };

	tn = omp_get_wtime();
	for (double t = t0; t < tMax; t += tau)
	{
		for (int i = 0; i < n; i++)
		{
			b[i] = -f(y, t, i);
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				memcpy(yy, y, sizeof(y));
				yy[j] += deltah;
				a[i][j] = (f(yy, t, i) - f(y, t, i)) / deltah;
			}
		}

		a[0][0] -= (double)(1.0 / tau);
		a[1][1] -= (double)(1.0 / tau);

		opred = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		if (!opred)
		{
			std::cout << "Определитель = 0!" << std::endl;
			exit(1);
		}

		opredn[0] = b[0] * a[1][1] - b[1] * a[0][1];
		opredn[1] = b[1] * a[0][0] - b[0] * a[1][0];

		for (int i = 0; i < n; i++)
		{
			p[i] = opredn[i] / opred;
		}

		for (int i = 0; i < n; i++)
		{
			y[i] += p[i];
		}
	}

	tk = omp_get_wtime();
	deltat = tk - tn;

	std::cout << "Неявный метод Эйлера" << std::endl;
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{
		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}

int main()
{
	EulerMethod(f, n, tMax, y);
	RK2(f, n, tMax, y);
	Correction(f, n, tMax, y);
	RK4(f, n, tMax, y);
	ImplictEulerMethod(f, n, tMax, y);
}