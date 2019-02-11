#ifndef DECK_PARAMETERS_H_
#define DECK_PARAMETERS_H_

#include <complex>

/******************************************************************************
 Namespace :: deck_parameters

 This module contains all parameters relevant for computing            
 the Ascoli Deck model.

 pi^-(p_a) + N(p_b,mu) --> pi^-(p_1) + p^-(p_2) + p^+(p_3) + N (p_d,mu')

 *****************************************************************************/

namespace Deck_parameters
{
  const double x_PI    = 3.1415926535897932384626433832795028841972;
  const double x_MPI   = 0.13957018;
  const double x_MPI0  = 0.1349766;
  const double x_MN    = 0.93827203;
  const double x_MK    = 0.493677;
  const double x_MK0   = 0.497614;
  const double x_HBC   = 0.19732;
  const std::complex<double> xi(0.0,1.0);
  const std::complex<double> xr(1.0,0.0);
  
}

#endif
