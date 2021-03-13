#include "../main/ReggeWheeler.h"

/*////////////////////////////////////
Circular Orbits in Schwarzschild
////////////////////////////////////*/

// angluar frequency for circular orbit
double Omega(double r, double M)
{
	double r_cube = r * r * r;
	return sqrt(M / r_cube);
}

// specific angular momentum of circular geodesic
double AngularMomentum(double r, double M)
{
	double r_square = r * r;
	return sqrt(r_square / (r - 3. * M));
}

// specific energy of circular geodesic
double Energy(double r, double M)
{
	return ((r - 2. * M) / sqrt(r * (r - 3. * M)));
}

// t-component of circular geodesic four-velocity, up and down
double ut_Down(double r, double M)
{
	return -Energy(r, M);
}

double ut_Up(double r, double M)
{
	return Energy(r, M) / f(r, M);
}

// phi-component of circular geodesic four-velocity, up and down
double uphi_Down(double r, double M)
{
	return AngularMomentum(r, M);
}

double uphi_Up(double r, double M)
{
	return utUp(r, M) * Omega(r, M);
}