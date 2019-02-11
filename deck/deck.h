#ifndef DECK_H_
#define DECK_H_

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <gsl/gsl_integration.h>

#include "deck_kinematics.h"
#include "rotational_functions.h"
#include "deck_pi_n.h"

typedef std::complex<double> cd;
typedef std::vector<double> vd;

namespace rf = rotational_functions;
namespace dk = deck_kinematics;
namespace dpN = deck_pi_n; 

namespace deck
{

  const double pi       = 3.1415926535897932384626433832795028841972;
  const cd xi(0.0,1.0);
  
  const double m_pi     = 0.13957018;
  const double m_p      = 0.93827203;
  const double m_rho    = 0.77;

  const double sq_m_pi  = m_pi * m_pi;
  const double sq_m_p   = m_p * m_p;
  const double sq_m_rho = m_rho * m_rho;

  const double E_beam   = 100.0;//190.0;

  const double s0       = sq_m_pi + sq_m_p + 2.0 * m_p * E_beam;

  const double t        = -0.1;


  struct arg_params {
    int J;
    int M;
    int S;
    int lam;
    double s;
    double sigma;
  };

  struct arg_piN_params{
    int M;
    double s;
    double sigma;
    double cos_th;
    int mu;
    int mu_p;
  };

  double lambda ( double x,
		  double y,
		  double z );

  double t_1 ( double s,
	       double sigma,
	       double cos_th );

  double cos_gamma ( double s,
		     double sigma );

  double cos_psi ( double s,
		   double sigma,
		   double t1 );

  double T_pi_p_0 ( double s,
		    double sigma,
		    double cos_th);


  double T_pi_p_1 ( double s,
		    double sigma,
		    double cos_th );

  double T_pi_p ( int M,
		  double s,
		  double sigma,
		  double cos_th );

  cd T_pi_p_dint ( int M,
		   double s,
		   double sigma,
		   double cos_th,
		   int mu,
		   int mu_p );

  double integrand_piN_real ( double phi,
                              void *p );
  
  double integrand_piN_imag ( double phi,
			      void *p );

  cd simple_piN ( double s,
                  double sigma,
                  double cos_th,
                  double phi,
                  int mu,
                  int mu_p );


  double t_min_max ( int pm,
		     double s );

  cd deck_bulk ( int lambda,
		 int S,
		 double s,
		 double sigma,
		 double cos_th );

  cd deck_reduced ( int M,
                    int lambda,
                    int S,
                    double s,
                    double sigma,
                    double cos_th );


  double intensity ( int lambda,
		     int S,
		     double s,
		     double sigma,
		     double cos_th );

  
  cd projection ( int J,
		  int M,
		  int S,
		  int lam,
		  double s,
		  double sigma );


  double integrand_real ( double cos_th,
			  void *p );

  double integrand_imag ( double cos_th,
			  void *p );
  

  void plot_wigner_d ( double J,
		       double m_p,
		       double m );

  void plot_2d_intensity ( double M3pi_i,
                           double M3pi_f,
                           double cos_th_i,
                           double cos_th_f,
                           int nstep_x,
                           int nstep_y );

  void plot_1d_intensity ( double cos_th_i,
			   double cos_th_f,
			   int nstep );

  void plot_1d_projection ( int J,
                            int M,
                            int S,
                            int lam,
                            double M3pi_i,
                            double M3pi_f,
                            int nstep);


  /****************************************************************************

   ***************************************************************************/


  double lambda ( double x,
                  double y,
                  double z )
  {
    return x * x + y * y + z * z - 2.0 * ( x * y + y * z + z * x );
  }

  double t_1 ( double s,
               double sigma,
               double cos_th )
  {
    double tmp1 = sigma + sq_m_pi -
      ( s + sigma - sq_m_pi ) * ( s + sq_m_pi - t ) / ( 2.0 * s );
    double tmp2 = sqrt ( lambda ( s, sigma, sq_m_pi ) *
			 lambda ( s, sq_m_pi, t ) ) / ( 2.0 * s ) * cos_th;
    return tmp1 + tmp2;
  }

  double cos_gamma ( double s,
                     double sigma )
  {
    double tmp1 = ( 2.0 * s * ( s0 + t - s - sq_m_p ) - 
		    ( s + sq_m_pi - t ) * ( s0 - s - sq_m_p ) );
    double tmp2 = sqrt ( lambda ( s, sq_m_pi, t ) * lambda ( s0, s, sq_m_p ) );
    return tmp1 / tmp2;
  }

  double cos_psi ( double s,
                   double sigma,
                   double t1 )
  {
    //    std::cout << "hey " << s << " " << sigma << " " << t1 << std::endl;
    double tmp1 = ( ( s + sigma - sq_m_pi ) * ( t1 + sigma - sq_m_pi ) - 
		    2.0 * sigma * ( sigma + t - 2.0 * sq_m_pi ) );
    double tmp2 = sqrt ( lambda ( s, sigma, sq_m_pi ) *
			 lambda ( sigma, t1, sq_m_pi ) );
    //    std::cout << tmp1 << " " << tmp2 << std::endl;
    return tmp1 / tmp2;
  }

  double T_pi_p_0 ( double s,
                    double sigma,
                    double cos_th )
  {
    double tmp1 = sq_m_pi + sq_m_p + 
      ( s0 - s - sq_m_p ) * ( s + sq_m_pi - sigma ) / ( 2.0 * s );
    double tmp2 = - sqrt ( lambda ( s, sigma, sq_m_pi ) 
			   * lambda ( s0, s, sq_m_p ) ) / ( 2.0 * s );
    tmp2        = tmp2 * cos_th * cos_gamma ( s, sigma );
    return tmp1 + tmp2;
  }

  double T_pi_p_1 ( double s,
                    double sigma,
                    double cos_th )
  {
    double tmp1 = - sqrt ( lambda ( s, sigma, sq_m_pi ) * 
			   lambda ( s0, s, sq_m_p ) ) / ( 2.0 * s );
    double tmp2 = sqrt ( 1.0 - cos_th * cos_th ) * 
      sqrt ( 1.0 - cos_gamma ( s, sigma ) * cos_gamma ( s, sigma ) );
    
    return tmp1 * tmp2;
  }

  double T_pi_p( int M,
		 double s,
		 double sigma,
		 double cos_th )
  {
    if ( M == 0 )
      {
	return T_pi_p_0 ( s, sigma, cos_th );
      }
    else if ( M == 1 )
      {
	return T_pi_p_1 ( s, sigma, cos_th );
      }
    else if ( M == -1 )
      {
	return T_pi_p_1 ( s, sigma, cos_th );
      }
  }


  cd T_pi_p_dint ( int M,
		   double s,
		   double sigma,
		   double cos_th,
		   int mu,
		   int mu_p )
  {
     gsl_integration_workspace *work_ptr =
       gsl_integration_workspace_alloc (10000);

     double lower_limit =  0.0;   /* lower limit a */
     double upper_limit =  2.0 * pi;   /* upper limit b */
     double abs_error = 0.01;   /* to avoid round-off problems */
     double rel_error = 1.0e-6;   /* the result will usually be much better */
     double result_real;          /* the result from the integration */
     double result_imag;
     double error;                /* the estimated error from the integration */
     gsl_function Freal;
     gsl_function Fimag;

     // abs_error and rel_error were 1.0e-8

     struct arg_piN_params params;
     params.M      = M;
     params.s      = s;
     params.sigma  = sigma;
     params.cos_th = cos_th;
     params.mu     = mu;
     params.mu_p   = mu_p;

     Freal.function = &integrand_piN_real;
     Fimag.function = &integrand_piN_imag;
     Freal.params = &params;
     Fimag.params = &params;
     
     std::cout << " phi - int success " << std::endl;

     gsl_integration_qags (&Freal, lower_limit, upper_limit, abs_error,
			   rel_error, 5000, work_ptr, &result_real, &error);
     gsl_integration_qags (&Fimag, lower_limit, upper_limit, abs_error,
			   rel_error, 5000, work_ptr, &result_imag, &error);
     std::cout << " phi - int after success " << std::endl;
     return cd( result_real, result_imag );
  }

  
  double integrand_piN_real ( double phi,
			      void *p )
  {
    struct arg_piN_params *params =(struct arg_piN_params *) p;
    int M         = (params->M);
    double s      = (params->s);
    double sigma  = (params->sigma);
    double cos_th = (params->cos_th);
    int mu        = (params->mu);
    int mu_p      = (params->mu_p);


    cd igrand;
    double W = std::sqrt (s);
    double MM = (double) M;
    //    igrand = simple_piN (s,sigma,cos_th,phi,mu,mu_p) * std::exp(-xi *MM * phi);
    igrand = dpN::amp_pi_N ( W, s0, t, sigma, cos_th, phi, mu, mu_p )
      * std::exp(-xi *MM * phi);
    std::cout << "igrand " << igrand << std::endl;
    return std::real ( igrand / 2.0 / pi );
  }


  double integrand_piN_imag ( double phi,
                              void *p )
  {
    struct arg_piN_params *params =(struct arg_piN_params *) p;
    int M         = (params->M);
    double s      = (params->s);
    double sigma  = (params->sigma);
    double cos_th = (params->cos_th);
    int mu        = (params->mu);
    int mu_p      = (params->mu_p);


    cd igrand;
    double W = std::sqrt (s);
    double MM = (double) M;
    //    igrand = simple_piN (s,sigma,cos_th,phi,mu,mu_p) * std::exp(-xi *MM * phi);
    igrand = dpN::amp_pi_N ( W, s0, t, sigma, cos_th, phi, mu, mu_p )
      * std::exp(-xi *MM * phi);
    std::cout << " igrand imag " << igrand << std::endl;
    return std::imag ( igrand / 2.0 / pi );
  }

  cd Deck_full ( double s,
		 double t,
		 double W,
		 double s1,
		 double z,
		 double phi )
  {
    int mu = 1;
    int mu_p = 1;
    cd piN = dpN::amp_pi_N ( W, s, t, s1, z, phi, mu, mu_p );
    double xt1 = dk::t_R ( W, s1, t, z );
    cd regge = 1.0 / ( xt1 - sq_m_pi );
    return regge * piN;
  }

  cd simple_piN ( double s,
		  double sigma,
		  double cos_th,
		  double phi,
		  int mu,
		  int mu_p )
  {
    if ( mu_p != mu )
      {
	cd amp(0.0,0.0);
	return amp;
      }
    else if ( mu_p == mu )
      {
	cd amp;
	double alpha = 9.0;
	double W = std::sqrt(s);
	double s_piN = dk::s_pi_N (W,sigma,s0,t,cos_th,phi) ;
	amp = s_piN ;//* std::exp ( -alpha * t );
	return amp;
      }
  }





  double t_min_max ( int pm,
		     double s )
  {
    double dpm  = (double) pm;
    double tmp1 = sq_m_pi + s - 
      ( s0 + sq_m_pi - sq_m_p ) * ( s0 + s - sq_m_p ) / ( 2.0 * s0 );
    double tmp2 = sqrt ( lambda ( s0, s, sq_m_p ) * 
			 lambda ( s0, sq_m_pi, sq_m_p ) ) / ( 2.0 * s0 );
    return tmp1 + tmp2 * dpm;
  }
  

  cd deck_bulk ( int lam,
                 int S,
                 double s,
                 double sigma,
                 double cos_th )
  {
    double xS       = (double) S;
    double xlam     = (double) lam;
    //double xt1      = t_1 ( s, sigma, cos_th );
    //    double xcos_psi = cos_psi ( s, sigma, xt1 );
    double W        = std::sqrt(s);
    double xt1      = dk::t_R ( W, sigma, t, cos_th );
    double psi      = dk::psi ( W, sigma, t, cos_th );

      //    double psi      =acos ( xcos_psi );
    cd v        = sqrt ( (2.0 * xS + 1.0 ) ) * 
      //      rf::djmn ( xS, 0.0, xlam, xcos_psi ) / ( sq_m_pi - xt1 );
      rf::wigner_d_matrix ( xS, 0.0, xlam, psi ) / ( sq_m_pi - xt1 );

    double b        = 1.7;
    double ff       = exp ( b * xt1 );
    return v ;//* ff;

  }

  cd deck_reduced ( int M,
		    int lam,
		    int S,
		    double s,
		    double sigma,
		    double cos_th )
  {
    cd db       = deck_bulk ( lam, S, s, sigma, cos_th ) ;
    //    double Tnuc = T_pi_p ( M, s, sigma, cos_th );
    cd Tnuc = T_pi_p_dint ( M, s, sigma, cos_th, 1, 1);

    return db * Tnuc;
  }

  
  double intensity ( int lam,
                     int S,
                     double s,
                     double sigma,
                     double cos_th )
  {
    cd BSlambda  = deck_bulk ( lam, S, s, sigma, cos_th );
    double Tpip0 = T_pi_p_0 ( s, sigma, cos_th );
    double Tpip1 = T_pi_p_1 ( s, sigma, cos_th );
    double v     = norm ( BSlambda );
    //    v            = v * ( Tpip0 * Tpip0 + 2.0 * Tpip1 * Tpip1 );
    double rho   = sqrt ( lambda ( s, sigma, sq_m_pi ) ) / ( 8.0 * pi * s );
    return v * rho;

  }


  cd projection ( int J,
                  int M,
                  int S,
                  int lam,
                  double s,
                  double sigma )
  {
    
    gsl_integration_workspace *work_ptr =
      gsl_integration_workspace_alloc (10000);

    double lower_limit = -1.0;   /* lower limit a */
    double upper_limit =  1.0;   /* upper limit b */
    double abs_error = 0.01;   /* to avoid round-off problems */
    double rel_error = 1.0e-6;   /* the result will usually be much better */
    double result_real;          /* the result from the integration */
    double result_imag;
    double error;                /* the estimated error from the integration */

    gsl_function Freal;
    gsl_function Fimag;

    // abs_error and rel_error were 1.0e-8  

    struct arg_params params;
    params.J     = J;
    params.M     = M;
    params.S     = S;
    params.lam   = lam;
    params.s     = s;
    params.sigma = sigma;
    
    Freal.function = &integrand_real;
    Fimag.function = &integrand_imag;
    Freal.params = &params;
    Fimag.params = &params;

    std::cout << " theta - int success " << std::endl;

    gsl_integration_qags (&Freal, lower_limit, upper_limit, abs_error, 
			  rel_error, 5000, work_ptr, &result_real, &error);
    gsl_integration_qags (&Fimag, lower_limit, upper_limit, abs_error, 
			  rel_error, 5000, work_ptr, &result_imag, &error);

    std::cout << " theta - int after success " << std::endl;
    return cd( result_real, result_imag );

  }


  double integrand_real ( double cos_th,
			  void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    int J        = (params->J);
    int M        = (params->M);
    int S        = (params->S);
    int lam      = (params->lam);
    double s     = (params->s);
    double sigma = (params->sigma);

    double xJ    = (double) J;
    double xM    = (double) M;
    double xlam  = (double) lam;
    
    double theta = acos(cos_th);
    double wd = rf::wigner_d_matrix ( xJ, xM, xlam, theta );

    cd dr = deck_reduced ( M, lam, S, s, sigma, cos_th );
    //    double wd = rf::djmn ( xJ, xM, xlam, cos_th );
    cd igrand = std::sqrt( 2.0 * xJ + 1.0 ) * dr * wd / 2.0;
    return std::real ( igrand );
  }

  double integrand_imag ( double cos_th,
                          void *p )
  {
    struct arg_params *params = (struct arg_params *) p;
    int J        = (params->J);
    int M        = (params->M);
    int S        = (params->S);
    int lam      = (params->lam);
    double s     = (params->s);
    double sigma = (params->sigma);

    double xJ   = (double) J;
    double xM   = (double) M;
    double xlam = (double) lam;

    double theta = acos(cos_th);
    double wd = rf::wigner_d_matrix ( xJ, xM, xlam, theta );

    cd dr     = deck_reduced ( M, lam, S, s, sigma, cos_th );
    //    double wd = rf::djmn ( xJ, xM, xlam, cos_th );
    cd igrand = std::sqrt( 2.0 * xJ + 1.0 ) * dr * wd / 2.0;
    return std::imag ( igrand );
  }



  void plot_wigner_d ( double J,
		       double m_p,
		       double m
		       )
  {
    double theta_i = 0.0;
    double theta_f = pi;
    int nstep = 1000;
    double step = ( theta_f - theta_i ) / ( (double) nstep );
    double wd = 0.0;
    double theta = theta_i;
    double cos_theta = 0.0;
    for ( int i = 1; i <= nstep; i++ )
      {
	cos_theta = std::cos( theta );
	//	wd = rf::djmn ( J, m_p, m, cos_theta );
	std::cout << theta << " " << wd << std::endl;
	theta = theta + step;
      }
    return;
  }
  
  void plot_2d_intensity ( double M3pi_i,
			   double M3pi_f,
			   double cos_th_i,
			   double cos_th_f,
			   int nstep_x,
			   int nstep_y )
  {

    int lam = 0;
    int S = 0;
    double sigma = sq_m_rho;

    double M3pi = M3pi_i;
    double cos_th = cos_th_i;
    double step_x = ( M3pi_f - M3pi_i ) / ( (double) nstep_x );
    double step_y = ( cos_th_f - cos_th_i ) / ( (double) nstep_y );
    double Intens = 0.0;
    double s = 0.0;
    for ( int j = 1; j <= nstep_x; j++ )
      {
	for ( int k = 1; k <= nstep_y; k++ )
	  {
	    s = M3pi * M3pi;
	    Intens = intensity ( lam, S, s, sigma, cos_th );
	    std::cout << M3pi << " " << cos_th << " " << Intens * M3pi << std::endl;
	    cos_th = cos_th + step_y;
	  }
	M3pi = M3pi + step_x;
	cos_th = cos_th_i;
	std::cout << " " << std::endl;
      }
    return;
  }


  void plot_1d_intensity ( double cos_th_i,
                           double cos_th_f,
                           int nstep )
  {
    int lam = 0;
    int S = 0;
    double sigma = sq_m_rho;
    double s = 5.0 * 5.0;

    double cos_th = cos_th_i;
    double step = ( cos_th_f - cos_th_i ) / ( (double) nstep );
    double Intens;
    for ( int i = 1; i <= nstep; i++ )
      {
	Intens = intensity ( lam, S, s, sigma, cos_th );
	std::cout << cos_th << " " << Intens << std::endl;
	cos_th = cos_th + step;
      }
    return;
  }


  void plot_1d_projection ( int J,
			    int M,
			    int S,
			    int lam,
			    double M3pi_i,
			    double M3pi_f,
			    int nstep )
    {
      double M3pi = M3pi_i;
      double step = ( M3pi_f - M3pi_i ) / ( (double) nstep );
      double sigma = sq_m_rho;
      double s = 0.0;
      for ( int i = 1; i <= nstep; i++ )
	{
	  s = M3pi * M3pi;
	  std::cerr << M3pi << std::endl;
	  //	  cd proj = projection ( J, M, S, lam, s, sigma );
	  //	  std::cerr << M3pi << " " << abs (proj) << std::endl;
	  //	  std::cout << M3pi << " " << abs (proj) << std::endl;
	  M3pi = M3pi + step;
	}

    }


  //  s0, t, s, sigma, cosGJ, phiGJ
  void read_production_mc ( vd &s_vec,
			    vd &t_vec,
			    vd &Wsq_vec,
			    vd &s1_vec,
			    vd &z_vec,
			    vd &phi_vec )
  {
    std::string MCfile = "production_MC.txt";


    std::ifstream inputFile(MCfile);

    double s_num, t_num, Wsq_num, s1_num, z_num, phi_num; 
    if (inputFile)
      {
        while ( inputFile >> s_num >> t_num >> 
		Wsq_num >> s1_num >> z_num >> phi_num  )
          {
            s_vec.push_back(s_num);
	    t_vec.push_back(t_num);
	    Wsq_vec.push_back(Wsq_num);
	    s1_vec.push_back(s1_num);
	    z_vec.push_back(z_num);
	    phi_vec.push_back(phi_num);
          }
      }
  }

  void out_deck ( )
  {
    vd s_vec, t_vec, Wsq_vec, s1_vec, z_vec, phi_vec;
    read_production_mc ( s_vec, t_vec, Wsq_vec, s1_vec, z_vec, phi_vec );
    
    cd dd;
    for ( int i = 0; i < s_vec.size(); i++ )
      {
	double WW = std::sqrt(Wsq_vec[i]);
	dd = Deck_full ( s_vec[i], t_vec[i], WW, 
			 s1_vec[i], z_vec[i], phi_vec[i] );

	std::cout << s_vec[i] << " " << t_vec[i] << " " 
		  << Wsq_vec[i] << " " << s1_vec[i] << " "
		  << z_vec[i] << " " << phi_vec[i] << " "
		  << std::real (dd) << " " << std::imag (dd) << std::endl;
      }
    
  }

}

#endif
