#pragma once

std::complex<double> RecurrenceHorizon(int nIn, int sIn, int lIn, double w, double M, std::complex<double> Anm1, std::complex<double> Anm2);
std::complex<double> RecurrenceInfinity(int nIn, int sIn, int lIn, double w, double M, std::complex<double> Anm1, std::complex<double> Anm2);
bool Summation(std::complex<double>* J, std::complex<double>* dJ, double factor, double dfactor, int s, int l, double w, double M, std::complex<double>(*Recurrence)(int, int, int, double, double, std::complex<double>, std::complex<double>), int Nmax = 0, int N = 0);

void HorizonBC(std::complex<double>* psi_final, std::complex<double>* dpsi_final, double* rfinal, int s, int l, double w, double M);
void InfinityBC(std::complex<double>* psi_final, std::complex<double>* dpsi_final, double* rfinal, int s, int l, double w, double M, double r0);