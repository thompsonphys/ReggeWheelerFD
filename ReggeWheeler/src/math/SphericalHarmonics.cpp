#include "../main/ReggeWheeler.h"

class Harmonic
{
public:
	int s;		  //Spin weight (if applicable)
	int l;		  //Orbital index
	int m;		  //Azimuthal index
	double theta; //Polar angle
	double phi;	  //Azimuthal angle
	std::complex<double> value;
	std::complex<double> dTheta;
	std::complex<double> d2Theta;
	std::complex<double> dPhi;
	std::complex<double> d2Phi;
	Harmonic(int s_in, int l_in, int m_in, double theta_in, double phi_in, bool ComputeDerivatives)
	{
		//Constructor for Harmonic class
		s = s_in;
		l = l_in;
		m = m_in;
		theta = theta_in;
		phi = phi_in;
		value = std::complex<double>(0., 0.);
		dTheta = value;
		d2Theta = value;
		dPhi = value;
		d2Phi = value;

		eval(ComputeDerivatives);
	}
	void eval(bool ComputeDerivatives)
	{
		if (s == 0)
		{
			value = Ylm(l, m, theta, phi);
		}
		else
		{
			value = Slm(s, l, m, theta, phi);
		}
		if (ComputeDerivatives)
		{
			compute_derivatives();
		}
	}
	void compute_derivatives()
	{
	}
};

//define associated legendre polynomials for all m
double assoc_legendre(int l, int m, double theta)
{
	if (m < 0)
	{
		int mabs = abs(m);
		return pow(-1, mabs) * factorial(l - mabs) / factorial(l + mabs) * std::sph_legendre(abs(l), mabs, theta);
	}
	else
	{
		return std::sph_legendre(abs(l), m, theta);
	}
}

//spherical harmonics
std::complex<double> Ylm(int l, int m, double theta, double phi)
{
	using namespace std::complex_literals;
	double mdbl = double(m);
	return assoc_legendre(l, m, theta) * exp(1i * mdbl * phi);
}

std::complex<double> Ylm(int s, int l, int m, double theta, double phi)
{
	assert(s == 0);
	using namespace std::complex_literals;
	double mdbl = double(m);
	return assoc_legendre(l, m, theta) * exp(1i * mdbl * phi);
}

//theta derivatives of spherical harmonics
std::complex<double> Ylm_dTheta(int l, int m, double theta, double phi)
{
	if (l == 0)
	{
		return 0.;
	}
	else
	{
		std::complex<double> term1 = m / tan(theta) * Ylm(l, m, theta, phi);
		if (m == l)
		{
			return term1;
		}
		else
		{
			using namespace std::complex_literals;
			return term1 + exp(-1i * phi) * sqrt((l - m) * (1 + l + m)) * Ylm(l, m + 1, theta, phi);
		}
	}
}

std::complex<double> Ylm_d2Theta(int l, int m, double theta, double phi)
{
	if (l == 0)
	{
		return 0.;
	}
	else
	{
		std::complex<double> term1 = m * (m / pow(tan(theta), 2) - 1. / pow(sin(theta), 2)) * Ylm(l, m, theta, phi);
		if (m == l)
		{
			return term1;
		}
		else
		{
			using namespace std::complex_literals;
			std::complex<double> term2 = exp(-1i * phi) * sqrt((l - m) * (1 + l + m)) * (1. + 2. * m) / tan(theta) * Ylm(l, m + 1, theta, phi);
			if (m == l - 1)
			{
				return term1 + term2;
			}
			else
			{
				return term1 + term2 + exp(-2i * phi) * sqrt((l - m - 1) * (l - m) * (l + m + 1) * (l + m + 2)) * Ylm(l, m + 2, theta, phi);
			}
		}
	}
}

//phi derivatives of spherical harmonics
std::complex<double> Ylm_dPhi(int l, int m, double theta, double phi)
{
	if (l == 0)
	{
		return 0.;
	}
	else
	{
		if (m == 0)
		{
			return 0.;
		}
		else
		{
			using namespace std::complex_literals;
			return 1i * double(m) * Ylm(l, m, theta, phi);
		}
	}
}

std::complex<double> Ylm_d2Phi(int l, int m, double theta, double phi)
{
	if (l == 0)
	{
		return 0.;
	}
	else
	{
		if (m == 0)
		{
			return 0.;
		}
		else
		{
			using namespace std::complex_literals;
			double dblm = double(m);
			return -dblm * dblm * Ylm(l, m, theta, phi);
		}
	}
}

//spheroidal harmonics
std::complex<double> Slm(int s, int l, int m, double theta, double phi)
{
	if (s == 0)
	{
		return Ylm(l, m, theta, phi);
	}
	else
	{
		//Add Code
		return std::complex<double>(0., 0.);
	}
}