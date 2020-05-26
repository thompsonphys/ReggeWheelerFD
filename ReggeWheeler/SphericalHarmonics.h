#pragma once

double assoc_legendre(int l, int m, double theta);

std::complex<double> Ylm(int l, int m, double theta, double phi);
std::complex<double> Ylm(int s, int l, int m, double theta, double phi);

std::complex<double> Ylm_dTheta(int l, int m, double theta, double phi);
std::complex<double> Ylm_d2Theta(int l, int m, double theta, double phi);
std::complex<double> Ylm_dPhi(int l, int m, double theta, double phi);
std::complex<double> Ylm_d2Phi(int l, int m, double theta, double phi);

std::complex<double> Slm(int s, int l, int m, double theta, double phi);