#include "ReggeWheeler.h"

/*////////////////////////////////////
Circular Orbits in Schwarzschild
////////////////////////////////////*/

// angluar frequency for circular orbit
double Omega(double r, double M)
{
	return sqrt(M / pow(r, 3));
}

// specific angular momentum of circular geodesic
double AngularMomentum(double r, double M)
{
	return sqrt(pow(r, 2) / (r - 3. * M));
}

// specific energy of circular geodesic
double Energy(double r, double M)
{
	return ( ( r - 2. * M ) / sqrt( r * ( r - 3. * M ) ) );
}

// t-component of circular geodesic four-velocity, up and down
double utDown(double r, double M)
{
	return -Energy(r, M);
}

double utUp(double r, double M)
{
	return Energy(r, M) / f(r, M);
}

// phi-component of circular geodesic four-velocity, up and down
double upDown(double r, double M)
{
	return AngularMomentum(r, M);
}

double upUp(double r, double M)
{
	return utUp(r, M) * Omega(r, M);
}