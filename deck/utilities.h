#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <complex>
#include <algorithm>
#include <cmath>
#include <vector>

typedef std::complex<double> cd;
typedef std::vector<double> vd;


namespace utilities
{
  void convert_complex ( cd a,
			 double &mag,
			 double &phi )
			 
  {
    double pi = 3.1415926535;
    mag = std::sqrt( std::abs ( a ) );
    double si2 = std::imag ( a ) / mag;
    double co2 = std::real ( a ) / mag;
    
    if ( ( si2 > 0.0 ) && ( co2 > 0.0 ) )
      {
	phi = std::atan ( si2 / co2 );
      }
    if ( ( si2 > 0.0 ) && ( co2 < 0.0 ) )
      { 
	phi = pi - std::atan( - si2 / co2 );
      } 
    if ( ( si2 < 0.0 ) && ( co2 < 0.0 ) )
      { 
	phi = pi + std::atan( si2 / co2 );
      } 
    if ( ( si2 < 0.0 ) && ( co2 > 0.0 ) )
      { 
	phi = 2.0 * pi - std::atan( -si2 / co2 );
      } 
    phi = phi * 180 / pi;
    return;
    
  }



}


#endif
