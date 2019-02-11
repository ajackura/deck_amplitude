#ifndef DECK_AMPLITUDES_H_
#define DECK_AMPLITUDES_H_


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "rotational_functions.h"
#include "deck_parameters.h"
#include "deck_kinematics.h"
#include "deck_pion_pion.h"
#include "deck_pion_nucleon.h"
#include "pion_pion_scattering.h"
#include "pipi.hpp" 


typedef std::complex<double> cd;
typedef std::vector<double> vd;
typedef std::vector<cd> cvd;


namespace rf = rotational_functions;
namespace dk = deck_kinematics;
namespace dpp = deck_pion_pion;
namespace dpn = deck_pion_nucleon;
namespace pp = pion_pion_scattering;

namespace deck_amplitudes
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
  const double sq_m_p   = m_p * m_p;
  const double sq_m_rho = m_rho * m_rho;

  struct arg_params {
    double W;
    double s1;
    double s;
    double t;
    int J;
    int M;
    int L;
    int S;
    int mu;
    int mu_p;
  };


  cd full_deck ( double W,
		 double s1,
		 double s,
		 double t,
		 double z,
		 double phi,
		 double z1,
		 double phi1,
		 int mu,
		 int mu_p )
  {
    double psi = dk::psi ( W, s1, t, z );
    double th1 = std::acos ( z1 );

    cvd pipi_amp = { iz, iz, iz };
    cvd PS = { iz, iz, iz };
    cd pipi_A = iz;
    for ( int S = 0; S <= 2; S++ )
      {
	cd sum = 0.0;
	for ( int lam = -S; lam <= S; lam++ )
	  {
	    double rS   = (double) S;
	    double rlam = (double) lam;
	    double dpsi = rf::wigner_d_matrix ( rS, 0.0, rlam, psi );
	    cd DS       = rf::wigner_D_matrix ( rS, rlam, 0.0, phi1, th1, 0.0);
	    sum         = sum + dpsi * std::conj ( DS );
	  }
	PS [S]          = std::pow ( -1.0, S ) * sum;
	/*	if ( ( S == 0) || (S == 2) )
	  {
	    pipi pipi_obj (0);
	    pipi_A = pipi_obj.GKPRY_partial_wave (S, s1 );

	  }else
	  {
	    pipi pipi_obj (1);
	    pipi_A = pipi_obj.GKPRY_partial_wave ( S, s1 );
	    }*/
	pipi_A = pp::pipi_amp ( S, s1 );
	pipi_amp[S] = pipi_A * PS[S];
      }
    cd tmp = 2.0 * pipi_amp[0] / 3.0;
      //      - 3.0 * pipi_amp[1]
      //      + 10.0 * pipi_amp[2] / 3.0;
    cd pion_propagator = dpp::standard_pion_pion ( W, s1, s, t, z );
    //cd piN_amp         = dpn::amp_pi_N ( W, s, t, s1, z, phi, mu, mu_p );
    cd piN_amp         = dpn::simple_piN ( W, s1, s, t, z, phi, mu, mu_p );
    return tmp * pion_propagator * piN_amp;
  }

  cd Famp_integrand ( double W,
		      double s1,
		      double s,
		      double t,
		      double z,
		      double phi,
		      int J,
		      int M,
		      int L,
		      int S,
		      int mu,
		      int mu_p )
  {
    double theta = std::acos ( z );
    double rJ    = (double) J;
    double rM    = (double) M;
    double rL    = (double) L;
    double rS    = (double) S;
    double rmu   = (double) mu;
    double rmu_p = (double) mu_p;
    cd sum       = iz;
    cd tmp       = iz;
    
    double psi = dk::psi ( W, s1, t, z );

    for ( int lam = -S; lam <= S; lam++ )
      {
	double rlam = (double) lam;
	tmp = rf::clebsch_gordan ( rL, 0.0, rS, rlam, rJ, rlam );
	tmp = tmp * rf::wigner_d_matrix ( rJ, rM, rlam, theta );
	tmp = tmp * rf::wigner_d_matrix ( rS, 0.0, rlam, psi );
	sum = sum + tmp;
      }
    double factor = std::sqrt ( ( 2.0 * rL + 1.0 ) 
				* ( 2.0 * rS + 1.0 ) );
    factor = factor * std::pow ( -1.0, rS );
    //    cd pion_propagator = dpp::standard_pion_pion ( W, s1, s, t, z );
    cd pion_propagator = dpp::regge_pion_pion ( W, s1, s, t, z );
    //    cd piN_amp         = dpn::amp_pi_N ( W, s, t, s1, z, phi, mu, mu_p );
    cd piN_amp = dpn:: M_pp_p_mm ( W, s, t, s1, z, phi );

    
    //cd piN_amp = dpn::simple_piN ( W, s1, s, t, z, phi, mu, mu_p );
    sum = sum * pion_propagator * piN_amp * std::exp ( -ii * rM * phi );
    cd pipi_amp = pp::pipi_amp ( S, s1 );
    return factor * pipi_amp * sum;
  }


  double real_Famp_integrand ( double *x, 
			       size_t dim, 
			       void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s1 = (params->s1);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double   z = x[0];
    double phi = x[1];
    cd Famp    = Famp_integrand ( W, s1, s, t, z, phi, J, M, L, S, mu, mu_p );
    return std::real ( Famp );
  }

  double imag_Famp_integrand ( double *x,
                               size_t dim,
                               void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s1 = (params->s1);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double   z = x[0];
    double phi = x[1];
    cd Famp    = Famp_integrand ( W, s1, s, t, z, phi, J, M, L, S, mu, mu_p );
    return std::imag ( Famp );
  }


  cd Famp ( double W,
	    double s1,
	    double s,
	    double t,
	    int J,
	    int M,
	    int L,
	    int S,
	    int mu,
	    int mu_p )
  {
    double resultReal, errorReal, resultImag, errorImag; 

    double xl[2] = { -1.0, 0.0 };
    double xu[2] = {  1.0, 2.0 * pi };

    const gsl_rng_type *T;
    gsl_rng *r;

    struct arg_params params;
    params.W    = W;
    params.s1   = s1;
    params.s    = s;
    params.t    = t;
    params.J    = J;
    params.M    = M;
    params.L    = L;
    params.S    = S;
    params.mu   = mu;
    params.mu_p = mu_p;

    gsl_monte_function GReal = { &real_Famp_integrand, 2, &params };
    gsl_monte_function GImag = { &imag_Famp_integrand, 2, &params };

    size_t calls = 500000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc (2);
    
    gsl_monte_vegas_integrate (&GReal, xl, xu, 2, 10000, r, state, 
			       &resultReal, &errorReal);
    gsl_monte_vegas_integrate (&GImag, xl, xu, 2, 10000, r, state,
                               &resultImag, &errorImag);
    //display_results ("vegas warm-up", result, error);

    std::cerr << "converging... " << std::endl;

    //    do
    //      {
        gsl_monte_vegas_integrate (&GReal, xl, xu, 2, calls / 5, r, state,
                                   &resultReal, &errorReal);
	gsl_monte_vegas_integrate (&GImag, xl, xu, 2, calls / 5, r, state,
                                   &resultImag, &errorImag);
	//	std::cout
	  //	  << "result = " << resultReal
	  //	  << " sigma = " << errorReal
	  //	  << " chisq/dof = " << state->chisq << std::endl;
	//      }
    //    while (fabs (state->chisq - 1.0) > 0.5);

    //    display_results ("vegas final", result, error);

    gsl_monte_vegas_free (state);

    return resultReal + ii * resultImag;
  }





  cd Fhel_integrand ( double W,
                      double s1,
                      double s,
                      double t,
                      double z,
                      double phi,
                      int J,
                      int M,
                      int lam,
                      int S,
                      int mu,
                      int mu_p )
  {
    double theta = std::acos ( z );
    double rJ    = (double) J;
    double rM    = (double) M;
    double rlam  = (double) lam;
    double rS    = (double) S;
    double rmu   = (double) mu;
    double rmu_p = (double) mu_p;
    cd sum       = iz;
    cd tmp       = iz;

    double psi = dk::psi ( W, s1, t, z );

    tmp = rf::wigner_d_matrix ( rJ, rM, rlam, theta );
    tmp = tmp * rf::wigner_d_matrix ( rS, 0.0, rlam, psi );
    sum = tmp;
    double factor = std::sqrt ( ( 2.0 * rS + 1.0 ) ) ;
    cd pion_propagator = dpp::standard_pion_pion ( W, s1, s, t, z );
    cd piN_amp = dpn::simple_piN ( W, s1, s, t, z, phi, mu, mu_p );
    //cd piN_amp = dpn::amp_pi_N ( W, s, t, s1, z, phi, mu, mu_p );
    sum = sum * pion_propagator * piN_amp * std::exp ( -ii * rM * phi );
    cd pipi_amp = ir;
    factor = factor * std::sqrt( 2.0 * rJ + 1.0 ) / 2.0;
    return factor * pipi_amp * sum / 2.0 / pi ;
  }



  ////

  double real_Fhel_integrand ( double *x,
                               size_t dim,
                               void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s1 = (params->s1);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double   z = x[0];
    double phi = x[1];
    cd Famp    = Fhel_integrand ( W, s1, s, t, z, phi, J, M, L, S, mu, mu_p );
    return std::real ( Famp );
  }

  double imag_Fhel_integrand ( double *x,
                               size_t dim,
                               void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s1 = (params->s1);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double   z = x[0];
    double phi = x[1];
    cd Famp    = Fhel_integrand ( W, s1, s, t, z, phi, J, M, L, S, mu, mu_p );
    return std::imag ( Famp );
  }




  cd Fhel ( double W,
            double s1,
            double s,
            double t,
            int J,
            int M,
            int lam,
            int S,
            int mu,
            int mu_p )
  {
    double resultReal, errorReal, resultImag, errorImag;

    double xl[2] = { -1.0, 0.0 };
    double xu[2] = {  1.0, 2.0 * pi };

    const gsl_rng_type *T;
    gsl_rng *r;

    struct arg_params params;
    params.W    = W;
    params.s1   = s1;
    params.s    = s;
    params.t    = t;
    params.J    = J;
    params.M    = M;
    params.L    = lam;
    params.S    = S;
    params.mu   = mu;
    params.mu_p = mu_p;

    gsl_monte_function GReal = { &real_Fhel_integrand, 2, &params };
    gsl_monte_function GImag = { &imag_Fhel_integrand, 2, &params };

    size_t calls = 500000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&GReal, xl, xu, 2, 10000, r, state,
                               &resultReal, &errorReal);
    gsl_monte_vegas_integrate (&GImag, xl, xu, 2, 10000, r, state,
                               &resultImag, &errorImag);


    std::cerr << "converging... " << std::endl;

    gsl_monte_vegas_integrate (&GReal, xl, xu, 2, calls / 5, r, state,
			       &resultReal, &errorReal);
    gsl_monte_vegas_integrate (&GImag, xl, xu, 2, calls / 5, r, state,
			       &resultImag, &errorImag);

    gsl_monte_vegas_free (state);

    return resultReal + ii * resultImag;
  }




  void plot_1d_projection ( double s,
			    double t,
			    double s1,
			    int mu,
			    int mu_p,
			    int J,
                            int M,
                            int L,
                            int S,
                            double W_i,
                            double W_f,
                            int Npts )
  {
    double W = W_i;
    double step = ( W_f - W_i ) / ( (double) Npts );
    for ( int i = 1; i <= Npts; i++ )
      {
	std::cerr << W << std::endl;
	cd proj = Fhel ( W, s1, s, t, J, M, L, S, mu, mu_p );
	std::cerr << W << " " << abs (proj) << std::endl;          
	std::cout << W << " " << abs (proj) << std::endl;          
	W = W + step;
      }

  }




  void read_production_mc ( vd &s_vec,
                            vd &t_vec,
                            vd &Wsq_vec,
                            vd &s1_vec,
                            vd &z_vec,
                            vd &phi_vec,
			    vd &z1_vec,
			    vd &phi1_vec )
  {
    //    std::string MCfile = "production_MC_t0.1_v2.txt";
    std::string MCfile = "production_MC_t0.1_1e6.txt";

    std::ifstream inputFile(MCfile);

    // s, t, Wsq, s1, z, phi, z1, phi1
    double s_num, t_num, Wsq_num, s1_num, z_num, phi_num, z1_num, phi1_num;
    if (inputFile)
      {
        while ( inputFile >> s_num >> t_num >>
                Wsq_num >> s1_num >> z_num >> 
		phi_num >> z1_num >> phi1_num )
          {
            s_vec.push_back(s_num);
            t_vec.push_back(t_num);
            Wsq_vec.push_back(Wsq_num);
            s1_vec.push_back(s1_num);
            z_vec.push_back(z_num);
            phi_vec.push_back(phi_num);
	    z1_vec.push_back(z1_num);
	    phi1_vec.push_back(phi1_num);
          }
      }
  }

  void out_deck ( )
  {
    vd s_vec, t_vec, Wsq_vec, s1_vec, z_vec, phi_vec, z1_vec, phi1_vec;
    read_production_mc ( s_vec, t_vec, Wsq_vec, s1_vec, z_vec, 
			 phi_vec, z1_vec, phi1_vec );

    cd dd;
    int mu = 1;
    int mu_p = 1;
    for ( int i = 0; i < s_vec.size(); i++ )
      {
        double WW = std::sqrt(Wsq_vec[i]);
	dd = full_deck ( WW, s1_vec[i], s_vec[i], t_vec[i],
			 z_vec[i], phi_vec[i], z1_vec[i],
			 phi1_vec[i], mu, mu_p );

	std::cerr << std::real (dd) << " " << std::imag (dd) << std::endl;
	std::cout << s_vec[i] << " " << t_vec[i] << " "
		  << Wsq_vec[i] << " " << s1_vec[i] << " "
		  << z_vec[i] << " " << phi_vec[i] << " "
		  << z1_vec[i] << " " << phi1_vec[i] << " "
		  << std::real (dd) << " " << std::imag (dd) << std::endl;
      }

  }






}


#endif
