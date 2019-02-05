#ifndef DECK_PI_N_H_
#define DECK_PI_N_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include "deck_kinematics.h"
#include "pion_nucleon_scattering.h"

namespace dk = deck_kinematics;
namespace pn = pion_nucleon_scattering;

namespace deck_pi_n
{
  /****************************************************************************
   Global Variables
  ****************************************************************************/
  double mpi    = Deck_parameters::x_MPI;
  double mn     = Deck_parameters::x_MN;
  double sq_mpi = mpi * mpi;
  double sq_mn  = mn * mn;

  const double pi       = 3.1415926535;
  const cd xi(0.0,1.0);
  
  /****************************************************************************
   Function Prototypes
  ****************************************************************************/

  double X ( double W,
	     double s,
	     double t );
  
  double Y ( double W,
             double s,
             double t );

  double U ( double W,
             double s,
             double t );

  double V ( double W,
             double s,
             double t );

  cd M_pp_p_mm ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi );

  cd M_pm_m_mp ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi );

  cd M_pp_m_mm ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi );
  
  cd M_pm_p_mp ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi );

  cd amp_pi_N ( double W,
		double s,
		double t,
		double s1, 
		double z,
		double phi,
		int mu,
		int mu_p );

  double cos2_ep_p_xi ( double W,
			double s,
			double t );

  double cos2_ep_m_xi (double W,
		       double s,
		       double t );

  double sin2_ep_p_xi (double W,
		       double s,
		       double t );

  double sin2_ep_m_xi (double W,
		       double s,
		       double t );


  void plot_M_pp_p_mm (double s,
		       double t,
		       double s1,
		       double z,
		       double phi,
		       double Wi,
		       double Wf,
		       int Npts );


  /****************************************************************************
   Functions
  ****************************************************************************/

  double X ( double W,
             double s,
             double t )
  {
    double Ed   = dk::Ed ( W, s );
    double Eb   = dk::Eb ( W, s, t );
    double tmp1 = std::sqrt ( ( Eb + mn ) * ( Ed + mn ) );
    double tmp2 = std::sqrt ( ( Eb - mn ) * ( Ed - mn ) );
    return tmp1 + tmp2;
  }

  double Y ( double W,
             double s,
             double t )
  {
    double Ed   = dk::Ed ( W, s );
    double Eb   = dk::Eb ( W, s, t );
    double tmp1 = std::sqrt ( ( Eb + mn ) * ( Ed + mn ) );
    double tmp2 = std::sqrt ( (Eb - mn) * ( Ed - mn ) );
    return tmp1 - tmp2;
  }

  double U ( double W,
             double s,
             double t )
  {
    double Ed   = dk::Ed ( W, s );
    double Eb   = dk::Eb ( W, s, t );
    double tmp1 = std::sqrt ( ( Eb + mn ) * ( Ed - mn ) );
    double tmp2 = std::sqrt ( (Eb - mn) * ( Ed + mn ) );
    return tmp1 + tmp2;
  }

  double V ( double W,
             double s,
             double t )
  {
    double Ed   = dk::Ed ( W, s );
    double Eb   = dk::Eb ( W, s, t );
    double tmp1 = std::sqrt ( ( Eb + mn ) * ( Ed - mn ) );
    double tmp2 = std::sqrt ( (Eb - mn) * ( Ed + mn ) );
    return tmp1 - tmp2;
  }


  cd M_pp_p_mm ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi )
  {
    double xY = Y ( W, s, t );
    double xX = X ( W, s, t );
    double xU = U ( W, s, t );
    
    double s_pi_N = dk::s_pi_N ( W, s1, s, t, z, phi );

    std::cout << "hey "  << s_pi_N << std::endl;

    cd A_piN; // pi- p -> pi- p = A+ + A-
    cd B_piN;

    double tR = dk::t_R ( W, s1, t, z );

    std::cout << "you "  << tR << std::endl;

    cv2d ABC(2, cvd(3,0));
    pn::Full_ABC ( s_pi_N, t, tR, ABC );
    A_piN = ABC[0][0] + ABC[1][0];
    B_piN = ABC[0][1] + ABC[1][1];

    std::cout << "A "  << A_piN << std::endl;
    std::cout << "B "  << B_piN << std::endl;
    
    double cosemx = cos2_ep_m_xi ( W, s, t );
    double cosepx = cos2_ep_p_xi ( W, s, t );
    double sinepx = sin2_ep_p_xi ( W, s, t );
    double sin_z  = std::sqrt ( 1.0 - z * z );

    std::cout << "xxx "  << cosemx << " " << cosepx  << std::endl;

    double E1     = dk::E1 ( W, s1 );
    double p1     = dk::p1 ( W, s1 );

    double tmp    = cosepx * z + sinepx * sin_z * cos ( phi );
    cd tmp1   = ( xY * A_piN + xX * E1 * B_piN ) * cosemx;
    cd tmp2   = - xU * p1 * B_piN * tmp;

    std::cout << " tmp " << tmp1 << " " << tmp2 << std::endl;
    
    return tmp1 + tmp2;

  }

  cd M_pm_m_mp ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi )
  {

    double xY = Y ( W, s, t );
    double xX = X ( W, s, t );
    double xV = V ( W, s, t );

    double s_pi_N = dk::s_pi_N ( W, s1, s, t, z, phi );

    cd A_piN; // pi- p -> pi- p = A+ + A-
    cd B_piN;

    double tR = dk::t_R ( W, s1, t, z );

    cv2d ABC(2, cvd(3,0));
    pn::Full_ABC ( s_pi_N, t, tR, ABC );
    A_piN = ABC[0][0] + ABC[1][0];
    B_piN = ABC[0][1] + ABC[1][1];

    double cosepx = cos2_ep_p_xi ( W, s, t );
    double sinepx = sin2_ep_p_xi ( W, s, t );
    double sinemx = sin2_ep_m_xi ( W, s, t );
    double sin_z  = std::sqrt ( 1.0 - z * z );

    double E1     = dk::E1 ( W,s1 );
    double p1     = dk::p1 ( W, s1 );

    double tmp    = sinepx * z - cosepx * sin_z * cos ( phi );
    cd tmp1   = - ( xX * A_piN + xY * E1 * B_piN ) * sinemx;
    cd tmp2   = - xV * p1 * B_piN * tmp;

    return tmp1 + tmp2;
  }

  cd M_pp_m_mm ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi )
  {
    double xU = U ( W, s, t );

    double s_pi_N = dk::s_pi_N ( W, s1, s, t, z, phi );

    cd B_piN; // pi- p -> pi- p = A+ + A-

    double tR = dk::t_R ( W, s1, t, z );

    cv2d ABC(2, cvd(3,0));
    pn::Full_ABC ( s_pi_N, t, tR, ABC );
    B_piN = ABC[0][1] + ABC[1][1];

    double sinemx = sin2_ep_m_xi ( W, s, t );
    double sin_z  = std::sqrt ( 1.0 - z * z );

    double p1     = dk::p1 ( W, s1 );

    return xU * p1 * B_piN * sinemx * sin_z * sin ( phi );
  }

  cd M_pm_p_mp ( double W,
		 double s,
		 double t,
		 double s1,
		 double z,
		 double phi )
  {
    double xV = V ( W, s, t );

    double s_pi_N = dk::s_pi_N ( W, s1, s, t, z, phi );


    cd B_piN; // pi- p -> pi- p = A+ + A-

    double tR = dk::t_R ( W, s1, t, z );

    cv2d ABC(2, cvd(3,0));
    pn::Full_ABC ( s_pi_N, t, tR, ABC );
    B_piN = ABC[0][1] + ABC[1][1];

    double cosemx = cos2_ep_m_xi ( W, s, t );
    double sin_z  = std::sqrt ( 1.0 - z * z );

    double p1     = dk::p1 ( W, s1 );

    return xV * p1 * B_piN * cosemx * sin_z * sin ( phi );
  }


  cd amp_pi_N ( double W,
		double s,
		double t,
		double s1,
		double z,
		double phi,
		int mu,
		int mu_p )
  {
    if ( ( mu_p == 1 ) && ( mu == 1 ) )
      {
	cd A = M_pp_p_mm ( W, s, t, s1, z, phi );
	cd C = M_pp_m_mm ( W, s, t, s1, z, phi );
	return A + xi * C;
      }
    else if ( ( mu_p == -1 ) && ( mu == -1 ) )
      { 
	cd A = M_pp_p_mm ( W, s, t, s1, z, phi );
	cd C = M_pp_m_mm ( W, s, t, s1, z, phi );
        return A - xi * C;
      } 

    else if ( ( mu_p == 1 ) && ( mu == -1 ) )
      { 
	cd B = M_pm_m_mp ( W, s, t, s1, z, phi );
	cd D = M_pm_p_mp ( W, s, t, s1, z, phi );
	return B + xi * D;
      } 

    else if ( ( mu_p == -1 ) && ( mu == 1 ) )
      { 
	cd B = M_pm_m_mp ( W, s, t, s1, z, phi );
	cd D = M_pm_p_mp ( W, s, t, s1,z, phi );
	return -B + xi * D;
      } 
  }

  double cos2_ep_p_xi ( double W,
                        double s,
                        double t )
  {
    double cos_ep = dk::cos_eps ( W, s, t );
    double cos_xi = dk::cos_xi ( W, s, t );
    double tmp1   = std::sqrt ( ( 1.0 + cos_ep ) * ( 1.0 + cos_xi ) );
    double tmp2   = std::sqrt ( ( 1.0 - cos_ep ) * ( 1.0 - cos_xi ) );
    return ( tmp1 - tmp2 ) / 2.0;
  }

  double cos2_ep_m_xi (double W,
                       double s,
                       double t )
  {
    double cos_ep = dk::cos_eps ( W, s, t );
    double cos_xi = dk::cos_xi ( W, s, t );
    double tmp1   = std::sqrt ( ( 1.0 + cos_ep ) * ( 1.0 + cos_xi ) );
    double tmp2   = std::sqrt ( ( 1.0 - cos_ep ) * ( 1.0 - cos_xi ) );
    return ( tmp1 + tmp2 ) / 2.0;
  }

  double sin2_ep_p_xi (double W,
                       double s,
                       double t )
  {
    double cos_ep = dk::cos_eps ( W, s, t );
    double cos_xi = dk::cos_xi ( W, s, t );
    double tmp1   = std::sqrt ( ( 1.0 - cos_ep ) * ( 1.0 + cos_xi ) );
    double tmp2   = std::sqrt ( ( 1.0 + cos_ep ) * ( 1.0 - cos_xi ) );
    return ( tmp1 + tmp2 ) / 2.0;
  }

  double sin2_ep_m_xi (double W,
                       double s,
                       double t )
  {
    double cos_ep = dk::cos_eps ( W, s, t );
    double cos_xi = dk::cos_xi ( W, s, t );
    double tmp1   = std::sqrt ( ( 1.0 - cos_ep ) * ( 1.0 + cos_xi ) );
    double tmp2   = std::sqrt ( ( 1.0 + cos_ep ) * ( 1.0 - cos_xi ) );
    return ( tmp1 - tmp2 ) / 2.0;
  }

  void plot_M_pp_p_mm ( double s,
			double t,
			double s1,
			double z,
			double phi,
			double Wi,
			double Wf,
			int Npts )
  {
    double W = Wi;
    double step = ( ( Wf - Wi ) / ( (double) Npts ) );
    for ( int i = 1; i <= Npts; i++ )
      {
	cd amp = M_pp_p_mm ( W, s, t, s1, z, phi );
	std::cout << W << " " << real(amp) << " " << imag(amp) << std::endl;
	W = W + step;
      }
  }
}

#endif
