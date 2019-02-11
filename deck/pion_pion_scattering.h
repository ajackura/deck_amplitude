#ifndef PION_PION_SCATTERING_H_
#define PION_PION_SCATTERING_H_

#include<iostream>
#include<sstream>
#include<fstream>
#include<cstring>
#include<cmath>
#include<vector>
#include<complex>

#include "rotational_functions.h"

namespace rf = rotational_functions;

typedef std::complex<double> cd;
typedef std::vector<double> vd;
typedef std::vector<vd> v2d;
typedef std::vector<v2d> v3d;
typedef std::vector<cd> cvd;
typedef std::vector<cvd> cv2d;
typedef std::vector<cv2d> cv3d;

namespace pion_pion_scattering
{
  const double pi    = 3.1415926535897932384626433832795028841972;
  const double twopi = 2.0 * pi;
  const cd ii (0.0,1.0);
  const cd ir (1.0,0.0);
  const cd iz (0.0,0.0);
  const double m_pi     = 0.13957018;
  const double m_p      = 0.93827203;
  const double m_rho    = 0.77;
  const double m_K      = 0.493677;
  const double m_pi0    = 0.1349766;
  const double m_K0     = 0.497614;
  const double m_pi_a   = ( m_pi + m_pi0 ) / 2.0;
  const double m_K_a    = ( m_K + m_K0 ) / 2.0;

  const double sq_m_pi   = m_pi * m_pi;
  const double sq_m_p    = m_p * m_p;
  const double sq_m_rho  = m_rho * m_rho;
  const double sq_m_K    = m_K * m_K;
  const double sq_m_K0   = m_K0 * m_K0;
  const double sq_m_pi0  = m_pi0 * m_pi0;
  const double sq_m_pi_a = m_pi_a * m_pi_a;
  const double sq_m_K_a  = m_K_a * m_K_a;

  cd cot_delta_S0 ( double s )
  {
    double B0 =  7.26;
    double B1 = -25.3;
    double B2 = -33.1;
    double B3 = -26.6;
    double z0sq = sq_m_pi;
    double s0 = 4.0 * sq_m_K;
    
    double k = std::sqrt ( s / 4.0 - sq_m_pi );
    cd w = ( std::sqrt ( s ) - std::sqrt ( (cd) s0 - s ) ) /
      ( std::sqrt ( s ) + std::sqrt ( (cd) s0 - s ) );
    std::cerr << " hello " << w << std::endl;
    cd tmp = sq_m_pi / ( s - z0sq / 2.0 );
    tmp = tmp * ( z0sq / ( m_pi * std::sqrt ( s ) ) + B0 + B1 * w 
		  + B2 * std::pow ( w, 2 ) + B3 * std::pow ( w, 3 ) ); 
    return std::sqrt ( s ) / 2.0 / k  * tmp;
    
  }

  double phase_space ( double s )
  {
    return std::sqrt ( 1.0 - 4.0 * sq_m_pi / s );
  }
  
  cd partial_wave_amplitude ( double s,
			      int L,
			      int I )
  {
    cd cot_delta;
    if ( L == 0 )
      {
	cot_delta =  cot_delta_S0 ( s );
      }
    else
      {
	cot_delta = 0.0;
      }
    double ps = phase_space ( s ) / 16.0 / pi ;
    cd amp = 1.0 / ( cot_delta - ii );
    return amp / ps;
  }

  double lambda ( double x,
		  double y,
		  double z )
  {
    return x * x + y * y + z * z - 2.0 * ( x * y + y * z + z * x );
  }


  cv2d pipi_S_mat ( double s )
  {
    //    cv2d amp(2, cvd(2) );
    const vd f   = {  0.1968, -0.0154 };
    const vd sp  = { -0.0074,  0.9828 }; 
    const v2d c0 = { {  0.0337, -0.2826 } , { -0.2826,  0.3010 } };
    const v2d c1 = { { -0.3185,  0.0918 } , {  0.0918, -0.5140 } };
    const v2d c2 = { { -0.0942,  0.1669 } , {  0.1669,  0.1176 } };
    const v2d c3 = { { -0.5927, -0.2082 } , { -0.2082,  0.5204 } };
    const v2d c4 = { {  0.1957, -0.1386 } , { -0.1386, -0.3977 } };
    const v2d a  = { {  0.1131,  0.0150 } , {  0.0150, -0.3210 } };
    v3d c( 5, v2d( 2, vd(2) ) );
    cv2d rhoM( 2, cvd(2) );
    cv2d M( 2, cvd(2) );
    cd rho1, rho2, detT, scale, tmp1, tmp2, tmp3;
    cv2d Tinv( 2, cvd(2) );
    cv2d Tamp( 2, cvd(2) );
    
    for ( int i = 0; i <= 1; i++ )
      {
	for ( int j = 0; j <= 1; j++ )
	  {
	    c[0][i][j] = c0[i][j];
	    c[1][i][j] = c1[i][j];
	    c[2][i][j] = c2[i][j];
	    c[3][i][j] = c3[i][j];
	    c[4][i][j] = c4[i][j];
	    rhoM[i][j] = iz;
	  }
      }
    rho1 = std::sqrt ( ir - 4.0 * sq_m_pi_a / s );
    rho2 = std::sqrt ( ir - 4.0 * sq_m_K_a / s );
    rhoM[0][0] = rho1;
    rhoM[1][1] = rho2;
    
    scale = s / 4.0 / sq_m_K_a - 1.0;
    
    tmp1 = iz;
    tmp2 = iz;
    tmp2 = iz;

    for ( int i = 0; i <= 1; i++ )
      {
	for ( int j = 0; j <= 1; j++ )
	  { 
	    tmp1 = a[i][j] / ( s - sp[0] + ii * 1e-10 );
	    tmp2 = f[i] * f[j] / ( s - sp[1] + ii * 1e-10 );
	    tmp3 = iz;
	    for ( int k = 0; k <= 4; k++ )
	      {
		tmp3 = tmp3 + c[k][i][j] * std::pow ( scale, k );
	      }
	    M[i][j] = tmp1 + tmp2 + tmp3;
	    Tinv[i][j] = ( M[i][j] - ii * rhoM[i][j] ) / 16.0 / pi;
	  } 
      }
    
    detT = Tinv[0][0] * Tinv[1][1] - Tinv[0][1] * Tinv[1][0];

    Tamp[0][0] =   Tinv[1][1] / detT;
    Tamp[0][1] = - Tinv[0][1] / detT;
    Tamp[1][0] = - Tinv[1][0] / detT;
    Tamp[1][1] =   Tinv[0][0] / detT;
   
    return Tamp;


  }



  cd pipi_S ( double s )
  {
    double a11 = 0.1131;
    vd c11 = { 0.0337, -0.3185, -0.0942, -0.5927 };
    double s0 = -0.0074;
    double f1_1 = 0.1968;
    double s1   = 0.9828;
    double scale = s / ( std::pow ( m_K + m_K0, 2 ) ) - 1.0;
    cd tmp1 = std::sqrt ( (cd) lambda ( s, sq_m_pi, sq_m_pi ) );
    cd tmp2 = std::sqrt ( (cd) lambda ( s, sq_m_pi0, sq_m_pi0 ) );
    cd rho11 = ( tmp1 + tmp2 ) / 2.0 / s;
    double tmp = 0.0;
    for ( int i = 0; i <=3; i++ )
      {
	tmp = tmp + c11 [i] * std::pow ( scale, i );
      }
    cd M11 = a11 / ( s - s0 ) + f1_1 * f1_1 / ( s1 - s ) + tmp;
    //std::cerr << " het " << a11 << " " << tmp << " " << rho11 << std::endl;
    cd amp = 1.0 * pi / ( M11 - ii * rho11 ); // 16.0
    return amp;
  }


  cd pipi_P ( double s )
  {
    double Mrho = 0.7755;
    double Gamrho = 0.149;
    double R = 5.0;
    double p, p0, Gam, g;
    p = std::sqrt ( lambda (s, sq_m_pi, sq_m_pi) ) / 2.0 / std::sqrt ( s );
    p0 = std::sqrt ( lambda (Mrho * Mrho, sq_m_pi, sq_m_pi) ) / 2.0 / Mrho;
    Gam = Gamrho * ( p * p / p0 / p0 ) * ( 1.0 + R * R * p0 * p0 ) 
      / ( 1.0 + R * R * p * p );
    Gam = Gam * ( p / p0 ) * ( Mrho / std::sqrt ( s ) );
    g = Mrho * std::sqrt ( 8.0 * pi * Gamrho / p0 );
    return g * g / ( Mrho * Mrho - s - ii * Mrho * Gam ); 
  }

  cd pipi_D ( double s )
  {
    double Mf2 = 1.2755;
    double Gamf2 = 0.1867;
    double R = 5.0;
    double p, p0, Gam, g;
    p = std::sqrt ( lambda (s, sq_m_pi, sq_m_pi) ) / 2.0 / std::sqrt ( s );
    p0 = std::sqrt ( lambda (Mf2 * Mf2, sq_m_pi, sq_m_pi) ) / 2.0 / Mf2;
    Gam = Gamf2 * std::pow ( ( p * p / p0 / p0 ) * ( 1.0 + R * R * p0 * p0 )
			     / ( 1.0 + R * R * p * p ), 2 );
    Gam = Gam * ( p / p0 ) * ( Mf2 / std::sqrt ( s ) );
    g = Mf2 * std::sqrt ( 8.0 * pi * Gamf2 / p0 );
    return g * g / ( Mf2 * Mf2 - s - ii * Mf2 * Gam );
  }

  
  cd pipi_amp ( int L,
		double s )
  {
    cd amp;
    cv2d Samp = pipi_S_mat ( s );
    switch ( L )
      {
      case 0: amp = Samp[0][0]; //pipi_S ( s );
	break;
      case 1: amp = pipi_P ( s );
	break;
      case 2: amp = pipi_D ( s );
	break;
      default: amp = iz;
	break;
      }
    return amp;
  }

  void plot_pipi_amp ( int L,
		       int I,
		       double mi,
		       double mf,
		       int Npts )
  {
    double m = mi;
    double step = ( mf - mi ) / ( (double) Npts );
    for ( int i = 1; i <= Npts; i++ )
      {
	double s = m * m;
	//cd amp = partial_wave_amplitude ( s, L , I );
	//	cd amp = pipi_S ( s );
	cd amp = pipi_amp ( L, s );
	std::cout << m << " " << std::real ( amp ) 
		  << " " << std::imag ( amp ) << std::endl;
	m = m + step;
      }
  }
  

}

#endif
 
