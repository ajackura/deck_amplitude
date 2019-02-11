// Self-contained implementation of the GKPY paramterization for low-energy pion scattering.
// Based on: 10.1103/PhysRevD.83.074004
//
// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PIPI_
#define _PIPI_

#include <complex>
#include <iostream>
#include <stdlib.h>
#include <cmath>

//------------------------------------------------------------------------------
// Pion Scattering Amplitude Object
class pipi
{
//-----------------------------------------------------------------------------
protected:
int qn_I;
//-----------------------------------------------------------------------------
public:
static const double conv; // degrees to radians conversion
static const std::complex<double> xr, xi; // unit imaginary and real

static const double mPi, mK, mEta, mRho, mF2; //masses
static const double sthPi, sthK, sthEta; // particle thresholds
//-----------------------------------------------------------------------------
pipi();
pipi(int i);
void print_I();

double conformal(double s, double s0);
double elastic_mom( double s, double sth);
double legendre(int l, double x);

double phase_shift(int l, double s);
double inelasticity(int l, double s);

std::complex<double> GKPRY_partial_wave(int l, double s);
std::complex<double> GKPRY_iso_amp(double s, double z);
};
//-----------------------------------------------------------------------------
#endif
