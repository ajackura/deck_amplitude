#ifndef DECK_KINEMATICS_H_
#define DECK_KINEMATICS_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include "deck_parameters.h"


/******************************************************************************
 Namespace :: deck_kinematics

 This module contains all kinematic functions relevant for computing
 the Ascoli Deck model.

 pi^-(p_a) + N(p_b,mu) --> pi^-(p_1) + p^-(p_2) + p^+(p_3) + N (p_d,mu')


 Dependencies --- deck_parameters.h

 
 *****************************************************************************/


namespace deck_kinematics
{
  /****************************************************************************
   Global Variables
  ****************************************************************************/
  double mpi    = Deck_parameters::x_MPI;
  double mn     = Deck_parameters::x_MN;
  double sq_mpi = mpi * mpi;
  double sq_mn  = mn * mn;

  /****************************************************************************
   Function Prototypes
  ****************************************************************************/

  double Ea ( double W,
	      double t );

  double Ed ( double W,
	      double s );

  double Eb ( double W,
	      double s,
	      double t );

  double E1 ( double W, 
	      double s1 );

  double pa ( double W,
	      double t );
  
  double pb ( double W,
	      double s,
	      double t );

  double pd ( double W,
	      double s );
  
  double p1 ( double W, 
	      double s1 );

  double cos_xi ( double W,
		  double s,
		  double t );

  double cos_eps ( double W,
		   double s, 
		   double t );

  double psi ( double W,
	       double s1,
	       double t,
	       double z );

  double t_R ( double W,
	       double s1,
	       double t,
	       double z );
  
  double s_pi_N ( double W,
		  double s1,
		  double s,
		  double t,
		  double z,
		  double phi );

  double q1 ( double s1 );

  double s2 ( double W,
	      double s1,
	      double z1 );

  double s3 ( double W,
	      double s1,
	      double z1 );
  
  double q2 ( double W,
	      double s1,
	      double z1 );

  double p2 ( double W,
	      double s1,
	      double z1 );

  double z2 ( double W,
	      double s1,
	      double z1 );

  double cos_theta12 ( double W,
		       double s1, 
		       double z1 );

  void plot_Eabd_Pabd ( double s,
			double t,
			double Wi,
			double Wf,
			int Npts );


  double Ea ( double s );

  double Eb ( double s );

  double Ed ( double s );

  double E1 ( double s, 
	      double s1 );

  double pa ( double s );
  
  double pb ( double s );

  double pd ( double s );
  
  double p1 ( double s, 
	      double s1 );

  double cos_alpha ( double s );

  double cos_beta ( double s );

  double psi ( double s,
	       double s1,
	       double z );

  double t1 ( double s,
	      double s1,
	      double z );
  
  double s4 ( double s,
	      double s1,
	      double z,
	      double phi );

  double q1 ( double s1 );

  double s2 ( double s,
	      double s1,
	      double z1 );

  double s3 ( double s,
	      double s1,
	      double z1 );
  
  double q2 ( double s,
	      double s1,
	      double z1 );

  double p2 ( double s,
	      double s1,
	      double z1 );

  double z2 ( double s,
	      double s1,
	      double z1 );

  double cos_theta12 ( double s,
		       double s1, 
		       double z1 );

  void plot_Eabd_Pabd ( double s,
                        double t,
			double Wi,
			double Wf,
                        int Npts );

  void plot_psi_tR ( double s,
                     double t,
                     double s1,
                     double Wi,
                     double Wf,
                     int Npts );

  /****************************************************************************
   Functions
  ***************************************************************************/
  double Ea ( double W,
	      double t )
  {
    return ( W * W + sq_mpi - t ) / 2.0 / W;
  }

  double Ed ( double W,
	      double s )
  {
    return ( s - W * W - sq_mn ) / 2.0 / W;
  }

  double Eb ( double W,
	      double s,
	      double t )
  {
    return Ed ( W, s ) + W - Ea ( W, t );
  }

  double E1 ( double W, 
	      double s1 )
  {
    return ( W * W + sq_mpi - s1 ) / 2.0 / W;
  }

  double pa ( double W,
	      double t )
  {
    return sqrt ( Ea ( W, t ) * Ea ( W, t ) - sq_mpi );
  }
  
  double pb ( double W,
	      double s,
	      double t )
  {
    return sqrt ( Eb ( W, s, t ) * Eb ( W, s, t ) - sq_mn );
  }

  double pd ( double W,
	      double s )
  {
    return sqrt ( Ed ( W, s ) * Ed ( W, s ) - sq_mn );
  }
  
  double p1 ( double W, 
	      double s1 )
  {
    double tmp1 = std::pow ( W - mpi, 2 ) - s1;
    double tmp2 = std::pow ( W + mpi, 2 ) - s1;
    return std::sqrt ( tmp1 * tmp2 ) / 2.0 / W;
  }

  double cos_xi ( double W,
		  double s,
		  double t )
  {
    double tmp = pd ( W, s ) * cos_eps ( W, s, t ) + pa ( W, t );
    return tmp / pb ( W, s, t );
  }

  double cos_eps ( double W,
		   double s, 
		   double t )
  {
    double sq_pb = pb ( W, s, t ) * pb ( W, s, t );
    double sq_pd = pd ( W, s ) * pd ( W, s );
    double sq_pa = pa ( W, t ) * pa ( W, t );
    return ( sq_pb - sq_pd - sq_pa ) / 2.0 / pd ( W, s ) / pa ( W, t );
  }

  double psi ( double W,
	       double s1,
	       double t,
	       double z )
  {
    double sin_z = std::sqrt ( 1.0 - z * z );
    double x = ( W - E1 ( W, s1 ) ) * pa ( W, t ) * z	\
               - p1 ( W, s1 ) * Ea ( W, t );
    double y = sqrt ( s1 ) * pa ( W, t ) * sin_z;
    return atan2 ( y, x );
  }

  double t_R ( double W,
	       double s1,
	       double t,
	       double z )
  {
    double tmp1 = 2.0 * pa ( W, t ) * p1 ( W, s1 ) * z;
    double tmp2 = sq_mpi + s1 - 2.0 * Ea ( W, t ) * ( W - E1 ( W, s1 ) );
    return tmp1 + tmp2;
  }
  
  double s_pi_N ( double W,
		  double s1,
		  double s,
		  double t,
		  double z,
		  double phi )
  {
    double sin_z   = sqrt ( 1.0 - z * z );
    double sin_eps = sqrt ( 1.0 - cos_eps ( W, s, t ) * cos_eps ( W, s, t ) );
    double cos_a   = cos_eps ( W, s, t ) * z + sin_eps * sin_z * cos ( phi );
    double tmp1    = 2.0 * pd ( W, s ) * p1 ( W, s1 ) * cos_a;
    double tmp2    = sq_mpi + sq_mn + 2.0 * Ed ( W, s ) * E1 ( W, s1 );
    return tmp2 - tmp1;
  }

  double q1 ( double s1 )
  {
    return sqrt ( s1 - 4.0 * sq_mpi ) / 2.0;
  }

  double s2 ( double W,
	      double s1,
	      double z1 )
  {
    double tmp = 2.0 * p1 ( W , s1 ) * q1 ( s1 ) * W * z1 / sqrt ( s1 );
    return ( W * W + 3.0 * sq_mpi - s1 ) / 2.0 + tmp;
  }

  double s3 ( double W,
	      double s1,
	      double z1 )
  {
    return W * W + 3.0 * sq_mpi - s1 - s2 ( W, s1, z1 );
  }
  
  double q2 ( double W,
	      double s1,
	      double z1 )
  {
    return ( s2 ( W, s1, z1 ) - 4.0 * sq_mpi ) / 2.0;
  }

  double p2 ( double W,
	      double s1,
	      double z1 )
  {
    double tmp1 = std::pow ( W - mpi, 2 ) - s2 ( W, s1, z1 );
    double tmp2 = std::pow ( W + mpi, 2 ) - s2 ( W, s1, z1 );
    return std::sqrt ( tmp1 * tmp2 ) / 2.0 / W;
  }

  double z2 ( double W,
	      double s1,
	      double z1 )
  {
    double tmp1 = W / std::sqrt ( s2 ( W, s1, z1 ) );
    double tmp2 = 4.0 * p2 ( W, s1, z1 ) * q2 ( W, s1, z1 ) * tmp1;
    return ( s1 - s3 ( W, s1, z1 ) ) / tmp2;
  }

  double cos_theta12 ( double W,
		       double s1, 
		       double z1 )
  {
    double tmp1 = std::sqrt ( 1.0 + std::pow ( p1 ( W, s1 ), 2 ) / s1 );
    double tmp2 = tmp1 * q1 ( s1 ) * z1 - p1 ( W, s1 ) / 2.0;
    return tmp2 / p2 ( W, s1, z1 );
  }


  void plot_Eabd_Pabd ( double s,
			double t,
			double Wi,
			double Wf,
			int Npts )
  {
    std::ofstream outputFile1("plot_Deck_Ea_s_" + std::to_string(s) + \
			      "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile2("plot_Deck_Eb_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile3("plot_Deck_Ed_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile4("plot_Deck_Pa_s_" + std::to_string(s) + \
			      "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile5("plot_Deck_Pb_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile6("plot_Deck_Pd_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile7("plot_Deck_cos_xi_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    std::ofstream outputFile8("plot_Deck_cos_ep_s_" + std::to_string(s) + \
                              "_t_" + std::to_string(t) + ".txt");
    double W    = Wi;
    double step = ( Wf - Wi ) / ( (double) Npts );
    for ( int n = 1; n < Npts; n++ )
      {
	outputFile1 << W << " " << Ea ( W, t ) << std::endl;
	outputFile2 << W << " " << Ed ( W, s ) << std::endl;
	outputFile3 << W << " " << Eb ( W, s, t ) << std::endl;
	outputFile4 << W << " " << pa ( W, t ) << std::endl;
	outputFile5 << W << " " << pd ( W, s ) << std::endl;
	outputFile6 << W << " " << pb ( W, s, t ) << std::endl;
	outputFile7 << W << " " << cos_xi ( W, s, t ) << std::endl;
	outputFile8 << W << " " << cos_eps ( W, s, t ) << std::endl;
	W  = W + step;
      }
    outputFile1.close();
    outputFile2.close();
    outputFile3.close();
    outputFile4.close();
    outputFile5.close();
    outputFile6.close();
    outputFile7.close();
    outputFile8.close();
    return;
  }


  void plot_psi_tR ( double s,
		     double t,
		     double s1,
		     double Wi,
		     double Wf,
		     int Npts )
  {
    
    double zi      = -1.0;
    double zf      =  1.0;
    double z       =  zi;
    double step_z  =  0.1;
    int nstep_z    =  (int) ( ( zf - zi ) / step_z );
    for ( int i = 0; i <= nstep_z; i ++ )
      {
	std::ofstream outputFile1("plot_Deck_psi_s_" + std::to_string(s) + \
				  "_t_" + std::to_string(t) + \
				  "_s1_" + std::to_string(s1) + \
				  "_z_" + std::to_string(z) + ".txt");
	std::ofstream outputFile2("plot_Deck_tR_s_" + std::to_string(s) + \
                                  "_t_" + std::to_string(t) + \
				  "_s1_" + std::to_string(s1) + \
                                  "_z_" + std::to_string(z) + ".txt");
	double W = Wi;
	double step = ( Wf - Wi ) / ( (int) Npts );
	for ( int n = 1; n <= Npts; n++ )
	  {
	    outputFile1 << W << " " << psi ( W, s1, t, z ) << std::endl;
	    outputFile2 << W << " " << t_R ( W, s1, t, z ) << std::endl;
	    W  = W + step;
	  }
	z = z + step_z;
      }
    
  }

}

#endif

