#ifndef DECK_OBSERVABLES_H_
#define DECK_OBSERVABLES_H_

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


namespace rf = rotational_functions;
namespace dk = deck_kinematics;
namespace dpp = deck_pion_pion;
namespace dpn = deck_pion_nucleon;
namespace pp = pion_pion_scattering;

namespace deck_observables
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


  cd direct_integrand ( double W,
			double s1,
			double s,
			double t,
			double za,
			double phia,
			double zb,
			double phib,
			int J,
			int M,
			int L,
			int S,
			int mu,
			int mu_p )
  {
    double th_a  = std::acos ( za );
    double th_b  = std::acos ( zb );
    double rJ    = (double) J;
    double rM    = (double) M;
    double rL    = (double) L;
    double rS    = (double) S;
    double rmu   = (double) mu;
    double rmu_p = (double) mu_p;
    cd suma      = iz;
    cd tmpa      = iz;
    cd sumb      = iz;
    cd tmpb      = iz;
    
    double psia = dk::psi ( W, s1, t, za );
    double psib = dk::psi ( W, s1, t, zb );

    for ( int lam = -S; lam <= S; lam++ )
      {
	double rlam = (double) lam;
	tmpa = rf::clebsch_gordan ( rL, 0.0, rS, rlam, rJ, rlam );
	tmpa = tmpa * rf::wigner_d_matrix ( rJ, rM, rlam, th_a );
	tmpa = tmpa * rf::wigner_d_matrix ( rS, 0.0, rlam, psia );
	suma = suma + tmpa;
	tmpb = rf::clebsch_gordan ( rL, 0.0, rS, rlam, rJ, rlam );
        tmpb = tmpb * rf::wigner_d_matrix ( rJ, rM, rlam, th_b );
        tmpb = tmpb * rf::wigner_d_matrix ( rS, 0.0, rlam, psib );
        sumb = sumb + tmpb;
      }
    double factor = std::sqrt ( ( 2.0 * rL + 1.0 ) 
				* ( 2.0 * rS + 1.0 ) );
    factor = factor * std::pow ( -1.0, rS );
    cd pion_propagatora = dpp::standard_pion_pion ( W, s1, s, t, za );
    cd piN_ampa = dpn::simple_piN ( W, s1, s, t, za, phia, mu, mu_p );
    suma = suma * pion_propagatora * piN_ampa * std::exp ( -ii * rM * phia );
    cd pion_propagatorb = dpp::standard_pion_pion ( W, s1, s, t, zb );
    cd piN_ampb = dpn::simple_piN ( W, s1, s, t, zb, phib, mu, mu_p );
    sumb = sumb * pion_propagatorb * piN_ampb * std::exp ( -ii * rM * phib );

    /*    cd pipi_amp = ir;


    if ( ( S == 0) || (S == 2) )
      {
	pipi pipi_obj (0);
	pipi_amp = pipi_obj.GKPRY_partial_wave (S, s1 );

      }else
      {
	pipi pipi_obj (1);
        pipi_amp = pipi_obj.GKPRY_partial_wave ( S, s1 );
      }

    //    std::cerr << pipi_amp << std::endl;

      //    cd pipi_amp = ir;
      */
    double iso;
    switch (S)
      {
      case 0: iso =   1.0 /3.0;
	break;
      case 1: iso = - 1.0 / 2.0;
	break;
      case 2: iso =   1.0 / 6.0;
	break;
      default: iso = 0.0;
	break;
      }
    cd pipi_amp = iso * pp::pipi_amp ( S, s1 );
    double p_s  = 2.0 * dk::p1 ( W, s1 ) * dk::q1( s1 ) / std::sqrt ( s1 );
				
    return p_s * factor * factor * pipi_amp * std::conj(pipi_amp)
      * suma * std::conj(sumb);
  }


  double real_integrand ( double *x, 
			  size_t dim, 
			  void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double  s1  = x[0];
    double  za  = x[1];
    double phia = x[2];
    double  zb  = x[3];
    double phib = x[4];
    cd Famp     = direct_integrand ( W, s1, s, t, za, phia, zb, phib,
				   J, M, L, S, mu, mu_p );
    return std::real ( Famp );
  }

  double imag_integrand ( double *x,
			  size_t dim,
			  void *p )
  {
    struct arg_params *params =(struct arg_params *) p;
    double W  = (params->W);
    double s  = (params->s);
    double t  = (params->t);
    int J     = (params->J);
    int M     = (params->M);
    int L     = (params->L);
    int S     = (params->S);
    int mu    = (params->mu);
    int mu_p  = (params->mu_p);

    double  s1  = x[0];
    double  za  = x[1];
    double phia = x[2];
    double  zb  = x[3];
    double phib = x[4];
    cd Famp     = direct_integrand ( W, s1, s, t, za, phia, zb, phib,
                                   J, M, L, S, mu, mu_p );
    return std::imag ( Famp );
  }


  double intensity ( double W,
		     double s,
		     double t,
		     int J,
		     int M,
		     int L,
		     int S,
		     int mu,
		     int mu_p )
  {
    struct arg_params params;
    params.W    = W;
    params.s    = s;
    params.t    = t;
    params.J    = J;
    params.M    = M;
    params.L    = L;
    params.S    = S;
    params.mu   = mu;
    params.mu_p = mu_p;


    double resultReal, errorReal, resultImag, errorImag; 

    double xl[5] = { 4.0 * sq_m_pi, -1.0, 0.0, -1.0, 0.0 };
    double xu[5] = { std::pow( W - m_pi ,2), 1.0, 2.0 * pi, 1.0, 2.0 * pi };

    const gsl_rng_type *T;
    gsl_rng *r;


    gsl_monte_function GReal = { &real_integrand, 5, &params };
    // gsl_monte_function GImag = { &imag_integrand, 5, &params };

    size_t calls = 500000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc (5);
    
    gsl_monte_vegas_integrate (&GReal, xl, xu, 5, 50000, r, state, 
			       &resultReal, &errorReal);
    //    gsl_monte_vegas_integrate (&GImag, xl, xu, 5, 10000, r, state,
    //                               &resultImag, &errorImag);
    //display_results ("vegas warm-up", result, error);

    std::cerr << "converging... " << std::endl;

    //    do
    //      {
        gsl_monte_vegas_integrate (&GReal, xl, xu, 5, calls / 5, r, state,
                                   &resultReal, &errorReal);
	//	gsl_monte_vegas_integrate (&GImag, xl, xu, 5, calls / 5, r, state,
	//                                   &resultImag, &errorImag);
	std::cerr
	  << "result = " << resultReal
	  << " sigma = " << errorReal
	  << " chisq/dof = " << state->chisq << std::endl;
	//      }
    //    while (fabs (state->chisq - 1.0) > 0.5);

    //    display_results ("vegas final", result, error);

    gsl_monte_vegas_free (state);

    return resultReal;// + ii * resultImag;
  }


  void plot_1d_intensity  ( double s,
			    double t,
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
	double it = intensity ( W, s, t, J, M, L, S, mu, mu_p );
	//	it = pi / std::pow( 4.0 * pi, 8);
	//	it = it * std::pow ( 0.19732 / 16.0 / m_p, 2 ) / 4.0;
	std::cerr << W << " " << it << std::endl;          
	std::cout << W << " " << it << std::endl;          
	W = W + step;
      }

  }








}


#endif
