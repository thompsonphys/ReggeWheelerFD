#include "ReggeWheeler.h"

class PPSourceCircularOrbitSchw {
private:
	std::complex<double> EA_compute()
	{
		return -16. * PI * f(r0, M) * Energy(r0, M) / r0 / r0 * std::conj(Ylm(l, m, PI / 2., 0.));
	}
	std::complex<double> EB_compute()
	{
		if (l==0 || m==0)
		{
			return 0.;
		}
		else 
		{
			using namespace std::complex_literals;
			double dblm = double(m);
			double dbll = double(l);
			double ll = dbll * (dbll + 1.);
			double r03 = r0 * r0 * r0;
			return 16. * PI * 1i * dblm * f(r0, M) * AngularMomentum(r0, M) / ll / r03 * std::conj(Ylm(l, m, PI / 2., 0.));
		}
	}
	std::complex<double> EC_compute()
	{
		if (l == 0)
		{
			return 0.;
		}
		else
		{
			using namespace std::complex_literals;
			double dblm = double(m);
			double dbll = double(l);
			double ll = dbll * (dbll + 1.);
			double r03 = r0 * r0 * r0;
			return -16. * PI * f(r0, M) * AngularMomentum(r0, M) / ll / r03 * std::conj(Ylm_dTheta(l, m, PI / 2., 0.));
		}
	}
	std::complex<double> EE_compute()
	{
		return -8. * PI * Omega(r0, M) * AngularMomentum(r0, M) / r0 / r0 * std::conj(Ylm(l, m, PI / 2., 0.));
	}
	std::complex<double> EF_compute()
	{
		if (l < 2)
		{
			return 0.;
		}
		else
		{
			double dblm = double(m);
			double dbll = double(l);
			double ll = dbll * (dbll + 1.);
			double lfacterm = 1. / ((dbll - 1.) * dbll * (dbll + 1.) * (dbll + 2.));
			return -16. * PI * lfacterm * Omega(r0, M)* AngularMomentum(r0, M) / r0 / r0 * ( ll - 2. * dblm * dblm ) * std::conj(Ylm(l, m, PI / 2., 0.));
		}
	}
public:
	int l;
	int m;
	double r0;
	double M;
	std::complex<double> EA;
	std::complex<double> EB;
	std::complex<double> EC;
	std::complex<double> ED;
	std::complex<double> EE;
	std::complex<double> EF;
	std::complex<double> EG;
	std::complex<double> EH;
	std::complex<double> EJ;
	std::complex<double> EK;
	PPSourceCircularOrbitSchw(int l_in, int m_in, double r0_in, double M_in )
	{
		l = l_in;
		m = m_in;
		r0 = r0_in;
		M = M_in;
		EA = EA_compute();
		EB = EB_compute();
		EC = EC_compute();
		ED = 0.;
		EE = EE_compute();
		EF = EF_compute();
		EG = 0.;
		EH = 0.;
		EJ = 0.;
		EK = 0.;
	}
};