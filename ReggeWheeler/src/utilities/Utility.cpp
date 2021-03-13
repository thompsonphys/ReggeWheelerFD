#include "ReggeWheeler.h"

int factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double n(int l)
{
	double dbl = double(l);
	return ((dbl - 1.) * (dbl + 2.) / 2.);
}
