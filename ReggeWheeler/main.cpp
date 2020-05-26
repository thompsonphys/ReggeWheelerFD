#include "ReggeWheeler.h"

int main()
{
	double r0 = 10.;
	double M = 1.0;
	int l = 20;
	int s = 2;
	int m = 2;
	double w = double(m) * Omega(r0, M);

	double rfinal;

	std::complex<double> psi = 0.;
	std::complex<double> dpsi = 0.;

	std::cout << psi << std::endl;
	std::cout << dpsi << std::endl;

	HorizonBC(&psi, &dpsi, &rfinal, s, l, w, M);

	std::cout << psi << std::endl;
	std::cout << dpsi << std::endl;
	std::cout << rfinal << std::endl;

	InfinityBC(&psi, &dpsi, &rfinal, s, l, w, M, r0);

	std::cout << psi << std::endl;
	std::cout << dpsi << std::endl;
	std::cout << rfinal << std::endl;

	std::cin.get();

	return 0;
}