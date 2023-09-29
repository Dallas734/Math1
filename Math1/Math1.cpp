#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

double F1(double* y, double t) 
{
	return 2 * (y[0] - y[0] * y[1]);
}

double F2(double* y, double t)
{
	return -(y[1] - y[0] * y[1]);
}

const int n = 2;
double (*p[n])(double*, double) = { F1, F2 };
const double tMax = 10.0;
double y[n] = { 1.0, 3.0 };
const double tau = 0.001;
const double t0 = 0.0;

void EulerMethod(double (*f[])(double*, double), int n, double tMax, double* yc, double tau, double t0)
{
	setlocale(LC_ALL, "Russian");

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
			yy[i] = y[i] + tau * f[i](y, t);
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

void RK2(double (*f[])(double*, double), int n, double tMax, double yc[], double tau, double t0)
{
	setlocale(LC_ALL, "Russian");

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
			yy[i] = y[i] + 0.5 * tau * f[i](y, t);
		}

		for (int i = 0; i < n; i++)
		{
			ff[i] = f[i](yy, t + 0.5 * tau);
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

void Correction(double (*f[])(double*, double), int n, double tMax, double yc[], double tau, double t0)
{
	setlocale(LC_ALL, "Russian");

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
			yy[i] = y[i] + tau * f[i](y, t);
		}

		for (int i = 0; i < n; i++)
		{
			ff[i] = y[i] + tau * (f[i](yy, t + tau) + f[i](y, t)) / 2.0;
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

void RK4(double (*f[])(double*, double), int n, double tMax, double yc[], double tau, double t0)
{
	setlocale(LC_ALL, "Russian");

	const int m = 4;
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
			R[0][i] = tau * f[i](y, t);
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + 0.5 * R[0][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[1][i] = tau * f[i](yy, t + 0.5 * tau);
			
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + 0.5 * R[1][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[2][i] = tau * f[i](yy, t + 0.5 * tau);
		}

		for (int i = 0; i < n; i++)
		{
			yy[i] = y[i] + R[2][i];
		}

		for (int i = 0; i < n; i++)
		{
			R[3][i] = tau * f[i](yy, t + tau);
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

void ImplictEulerMethod(double (*f[])(double*, double), int n, double tMax, double yc[], double tau, double t0)
{
	setlocale(LC_ALL, "Russian");

	double t = t0;
	double* y = new double[n] {0.0};
	memcpy(y, yc, n * sizeof yc);
	double* yy = new double[n] {0.0};
	double tn, tk, deltat;
	double deltah = 0.0001, opred; 
	double* opredn = new double[n]{ 0.0 };
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n] {0.0};
	}
	double* b = new double[n] {0.0}; 
	double* p = new double[n] {0.0};

	tn = omp_get_wtime();
	for (double t = t0; t < tMax; t += tau)
	{
		for (int i = 0; i < n; i++)
		{
			b[i] = -f[i](y, t);
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				memcpy(yy, y, sizeof(y) * n);
				yy[j] += deltah;
				a[i][j] = (f[i](yy, t) + b[i]) / deltah;
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
	EulerMethod(p, n, tMax, y, tau, t0);
	RK2(p, n, tMax, y, tau, t0);
	Correction(p, n, tMax, y, tau, t0);
	RK4(p, n, tMax, y, tau, t0);
	ImplictEulerMethod(p, n, tMax, y, tau, t0);
}