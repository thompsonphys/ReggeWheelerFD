#include "ReggeWheeler.h"

double f(double r, double M)
{
	return (1.0 - 2.0 * M / r);
}

double rstar(double r, double M)
{
	return (r + 2. * M * log(r / (2. * M) - 1.));
}

double rFromRstar(double rstar, double M)
{
	double x, dx, xstar = rstar / (2. * M) - 1.;

	if (xstar > 4. / 10.)
	{
		x = xstar - log(xstar - log(xstar));
		int jmax = 12;
		int j = 0;
		double xlast = x;

		do
		{
			double dxstar = xstar - x - log(x);
			dx = dxstar * x * (1. + dxstar / (2. * pow((x + 1.), 2))) / (x + 1.);
			xlast = x;
			x += dx;
			j++;
		} while (xlast != x && j != jmax);

		return 2. * M * x + 2. * M;
	}
	else
	{
		double y = exp(xstar);
		x = y;
		int jmax = 13;
		int j = 0;
		double xlast = x;

		x = y;

		do
		{
			xlast = x;
			dx = (y * exp(-x) - x) / (x + 1.);
			x += dx;
			j++;
		} while ((xlast != x) && (j <= jmax));

		return 2. * M * x + 2. * M;
	}

}

double rMinus2M(double r, double M)
{
	return (r - 2. * M);
}

double ReggeWheelerPotential(int s, int l, double r, double M)
{
	double dbls = double(s);
	double dbll = double(l);

	return f(r, M) / r / r * (dbll * (dbll + 1.) + 2. * M * (1. - dbls * dbls) / r);
}