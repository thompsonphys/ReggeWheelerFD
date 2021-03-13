#pragma once

// angluar frequency for circular orbit
double Omega(double r, double M);

// specific angular momentum of circular geodesic
double AngularMomentum(double r, double M);

// specific energy of circular geodesic
double Energy(double r, double M);

// t-component of circular geodesic four-velocity, up and down
double ut_Down(double r, double M);

double ut_Up(double r, double M);

// phi-component of circular geodesic four-velocity, up and down
double uphi_Down(double r, double M);

double uphi_Up(double r, double M);