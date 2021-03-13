#include "../main/ReggeWheeler.h"

std::complex<double> RecurrenceHorizon(int nIn, int sIn, int lIn, double w, double M, std::complex<double> Anm1, std::complex<double> Anm2)
{
	using namespace std::complex_literals;

	double l = double(lIn);
	double s = double(sIn);
	double n = double(nIn);
	double l2 = l * l;
	double n2 = n * n;
	double s2 = s * s;

	double term1, term2;

	std::complex<double> denominator;

	term1 = 1. + l + l2 - 2. * n + 2. * n2 - s2;
	term2 = -1. + 2. * n - n2 + s2;

	denominator = n * (n - 4i * M * w);

	return (term1 * Anm1 + term2 * Anm2) / denominator;
}

std::complex<double> RecurrenceInfinity(int nIn, int sIn, int lIn, double w, double M, std::complex<double> Anm1, std::complex<double> Anm2)
{
	using namespace std::complex_literals;

	double l = double(lIn);
	double s = double(sIn);
	double n = double(nIn);
	double l2 = l * l;
	double n2 = n * n;
	double s2 = s * s;

	std::complex<double> term1, term2, denominator;

	term1 = 1i * (l + l2 + n - n2);
	term2 = 1i * 2. * M * w * (1. - 2. * n + n2 - s2);
	;

	denominator = 2. * n;

	return (term1 * Anm1 + term2 * Anm2) / denominator;
}

bool Summation(std::complex<double> *J, std::complex<double> *dJ, double factor, double dfactor, int s, int l, double w, double M, std::complex<double> (*Recurrence)(int, int, int, double, double, std::complex<double>, std::complex<double>), int Nmax, int N)
{
	int n;
	std::complex<double> An = 1., Anm1 = 0., Anm2 = 0., lastJ = 0., lastdJ = 0., increment = 0.;
	double lastIncrement = 1.E40;

	*J = 0.;
	*dJ = 0.;

	double fN = pow(factor, N);

	for (n = N; n <= N + Nmax; n++)
	{
		increment = An * fN;

		lastJ = *J;
		lastdJ = *dJ;

		*J += increment;
		*dJ += double(n) * dfactor * increment / factor;

		if (n > N + 4 && *J == lastJ && *dJ == lastdJ)
		{
			return true;
		}

		Anm2 = Anm1;
		Anm1 = An;
		An = (*Recurrence)(n + 1, s, l, w, M, Anm1, Anm2);
		fN *= factor;
	}

	return false;
}

void HorizonBC(std::complex<double> *psi_final, std::complex<double> *dpsi_final, double *rfinal, int s, int l, double w, double M)
{
	using namespace std::complex_literals;

	std::complex<double> J = 0., dJ = 0.;
	double factor = 0., dfactor = 0., r, r2, drstardr, step;

	int k = 0;
	int kmax = 3;

	// we initialize at r = 2M + 1.E-5
	step = 1.E-5;
	r = 2. * M + step;
	factor = f(r, M);
	r2 = r * r;
	dfactor = 2. * M / r2;

	while (!Summation(&J, &dJ, factor, dfactor, s, l, w, M, RecurrenceHorizon, 10, 0) && k < kmax)
	{
		step /= 1.5;
		r = 2. * M + step;
		factor = f(r, M);
		r2 = r * r;
		dfactor = 2. * M / r2;
		k++;
	}

	drstardr = 1. / f(r, M);

	std::complex<double> psi_exp = exp(1i * w * rstar(r, M));
	std::complex<double> dpsi_exp = 1i * w * drstardr * psi_exp;

	*rfinal = r;
	*psi_final = J * psi_exp;
	*dpsi_final = J * dpsi_exp + dJ * psi_exp;
}

void InfinityBC(std::complex<double> *psi_final, std::complex<double> *dpsi_final, double *rfinal, int s, int l, double w, double M, double r0)
{
	using namespace std::complex_literals;

	std::complex<double> J = 0., dJ = 0.;
	double factor = 0., dfactor = 0., reval, drstardr;

	int k = 0;
	int kmax = 10;

	// we initialize at r = r0
	reval = r0 + 2. * PI / w;
	dfactor = -w;

	do
	{
		reval += 2. * PI / w;
		factor = 1. / (w * reval);
		k++;
	} while (!Summation(&J, &dJ, factor, dfactor, s, l, w, M, RecurrenceInfinity, 75, 0) && k <= kmax);

	drstardr = 1. / f(reval, M);

	std::complex<double> psi_exp = exp(1i * w * rstar(reval, M));
	std::complex<double> dpsi_exp = 1i * w * drstardr * psi_exp;

	*rfinal = reval;
	*psi_final = J * psi_exp;
	*dpsi_final = J * dpsi_exp + dJ * psi_exp;
}
