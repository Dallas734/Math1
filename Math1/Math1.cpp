#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

double f(double* y, double, int i)
{
	double w;

	switch (i)
	{
	case 0: return w = 2 * (y[0] - y[0] * y[1]); break;
	case 1: return w = -(y[1] - y[0] * y[1]); break;
	}
}

void EulerMethod()
{
	const int n = 2;
	const double t0 = 0;
	const double tMax = 100.0;
	const double tau = 0.01;
	double t = t0;
	double y[n] = { 1.0, 3.0 };
	double yy[n] = { 0.0 };
	double tn, tk, deltat;

	setlocale(LC_ALL, "Russian");
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
	std::cout << "Прошедшее время: " << deltat << std::endl;

	for (int i = 0; i < n; ++i)
	{
		std::cout << "[" << i << "]: " << y[i] << std::endl;
	}
}



int main()
{
	
}