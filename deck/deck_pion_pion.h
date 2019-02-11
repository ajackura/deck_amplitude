#ifndef DECK_PION_PION_H_
#define DECK_PION_PION_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>

#include "deck_parameters.h"
#include "deck_kinematics.h"


typedef std::complex<double> cd;
typedef std::vector<double> vd;


namespace dk = deck_kinematics;

namespace deck_pion_pion
{

  const double pi    = 3.1415926535897932384626433832795028841972;
  const double twopi = 2.0 * pi;
  const cd ii (0.0,1.0); // complex 'i'
  const cd ir (1.0,0.0); // complex '1'
  const cd iz (0.0,0.0); // complex '0'
  const double m_pi     = 0.13957018;
  const double m_p      = 0.93827203;
  const double m_rho    = 0.77;

  const double sq_m_pi  = m_pi * m_pi;


  cd standard_pion_pion ( double W,
			  double s1,
			  double s,
			  double t,
			  double z )
  {
    double tR = dk::t_R ( W, s1, t, z );
    return ir * 1.0 / ( sq_m_pi - tR );
  }

  cd regge_pion_pion ( double W,
		       double s1,
		       double s,
		       double t,
		       double z )
  {
    double s0    = 1.0;
    double tR    = dk::t_R ( W, s1, t, z );
    double alpha = tR - sq_m_pi;
    double s_m_u = W * W - sq_m_pi + 0.5 * ( tR - t - s1 );
    cd regge     = std::exp ( -ii * pi * alpha / 2.0 );
    regge        = regge * std::pow ( s_m_u / 2.0 / s0, alpha );
    return regge / ( sq_m_pi - tR );
  }

}

#endif
